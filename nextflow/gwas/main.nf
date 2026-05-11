// FinnGen regenie GWAS pipeline — Nextflow DSL2 port of wdl/gwas/*.wdl.
//
// Three variants are covered by the same pipeline, selected via params files:
//   1) Base release GWAS over bgen     (params/r13_binary_release.json)
//   2) Chip GWAS over bgen             (params/r13_chip_binary.json)
//   3) Chip GWAS over pgen             (params/r13_chip_pgen.json)
//
// Feature flags that differentiate the three variants:
//   input_format            bgen | pgen
//   step2_resource_profile  dynamic | fixed
//   step2_strip_chr         false | true      (pgen: strip 'chr' from CHROM)
//   gather_sort_flag        '-k1,1g -k2,2g' | '-k1,1V -k2,2g'
//   skip_qqplot             false | true
//   skip_summary            false | true      (chip variants)
//   skip_coding_gather      false | true      (chip variants)
//   run_sex_specific        true  | false

nextflow.enable.dsl = 2

include { STEP1 }         from './modules/step1.nf'
include { STEP2 }         from './modules/step2.nf'
include { GATHER }        from './modules/gather.nf'
include { SUMMARY }       from './modules/summary.nf'
include { CODING_GATHER } from './modules/coding_gather.nf'

workflow {

    // --- Pheno chunks -----------------------------------------------------
    // phenolist is a TSV: each row = one pheno chunk, columns = pheno names.
    // Chunk id = 8-char md5 of the joined name list — stable across resumes.
    pheno_chunks_ch = Channel.fromPath(params.phenolist)
        .splitCsv(sep: '\t', strip: true)
        .map { row -> row.findAll { it } }
        .filter { it }
        .map { row ->
            def joined = row.join(',')
            def id = joined.md5().take(8)
            tuple(id, row)
        }

    // --- GRM trio ---------------------------------------------------------
    // Stage .bed/.bim/.fam as three distinct paths so Nextflow can localise
    // them together in the task workdir.
    grm_prefix = params.grm_bed - ~/\.bed$/
    grm_ch = Channel.value(
        tuple(
            file("${grm_prefix}.bed"),
            file("${grm_prefix}.bim"),
            file("${grm_prefix}.fam"),
        )
    )

    // --- Shared cov_pheno -------------------------------------------------
    cov_pheno_ch = Channel.value(file(params.cov_pheno))

    // --- STEP 1 -----------------------------------------------------------
    if (params.step1_results) {
        // Reuse published step1 outputs — skip running STEP1 entirely.
        def grm_base = file(params.grm_bed).getBaseName()
        step1_resolved = pheno_chunks_ch.map { chunk_id, phenolist ->
            def dir = file("${params.step1_results}/${chunk_id}")
            def pred       = file("${dir}/${grm_base}.${chunk_id}.pred.list")
            def loco       = phenolist.collect { pheno -> file("${dir}/${grm_base}.${pheno}.loco.gz") }
            def firth_list = file("${dir}/${grm_base}.${chunk_id}.firth.list")
            def firth_gz   = params.is_binary
                ? phenolist.collect { pheno -> file("${dir}/${grm_base}.${pheno}.firth.gz") }
                : []
            def covars_used = file("${dir}/new_covars").text.trim()
            def single_sex  = file("${dir}/is_single_sex").text.trim() == 'true'
            tuple(chunk_id, phenolist, pred, loco, firth_list, firth_gz, covars_used, single_sex)
        }
    } else {
        step1_input_ch = pheno_chunks_ch
            .combine(grm_ch)
            .combine(cov_pheno_ch)
            .map { chunk_id, phenolist, bed, bim, fam, cov_pheno ->
                tuple(chunk_id, phenolist, bed, bim, fam, cov_pheno)
            }

        STEP1(step1_input_ch)

        // Resolve the new_covars / is_single_sex sidecar files into values so
        // STEP2 can consume them without re-staging.
        step1_resolved = STEP1.out.main.map { row ->
            def (chunk_id, phenolist, pred, loco, firth_list, firth_gz, new_covars_file, is_single_sex_file, log_file) = row
            def covars_used = new_covars_file.text.trim()
            def single_sex  = is_single_sex_file.text.trim() == 'true'
            tuple(chunk_id, phenolist, pred, loco, firth_list, firth_gz, covars_used, single_sex)
        }
    }

    // --- Genotype chunk list ---------------------------------------------
    // bgen variant: <path>.bgen + <path>.bgen.bgi + <path>.bgen.sample siblings.
    // pgen variant: <prefix>.pgen + <prefix>.pvar + <prefix>.psam siblings.
    // Param name stays 'bgenlist' for backwards-compat with the r13 JSONs.
    geno_ch = Channel.fromPath(params.bgenlist)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { p ->
            if (params.input_format == 'pgen') {
                def prefix = p.replaceFirst(/\.pgen$/, '')
                tuple(file(p), file("${prefix}.pvar"), file("${prefix}.psam"))
            } else {
                tuple(file(p), file("${p}.bgi"), file("${p}.sample"))
            }
        }

    // --- STEP 2 -----------------------------------------------------------
    if (params.step2_results) {
        // Reuse published step2 outputs — skip running STEP2 entirely.
        all_regenie_ch = Channel.fromPath("${params.step2_results}/*.regenie.gz").toList()
        all_sex_ch     = Channel.fromPath("${params.step2_results}/*.sex_spec.gz").toList()
        phenos_flat    = pheno_chunks_ch.flatMap { chunk_id, phenolist -> phenolist }

        regenie_by_pheno = phenos_flat
            .combine(all_regenie_ch)
            .map { pheno, files ->
                def matched = files.findAll { it.getName().endsWith("_${pheno}.regenie.gz") }
                tuple(pheno, matched)
            }

        sex_by_pheno = phenos_flat
            .combine(all_sex_ch)
            .map { pheno, files ->
                def matched = files.findAll { f ->
                    def n = f.getName()
                    n.endsWith(".${pheno}.sex_spec.gz") || n.contains('NOT_DONE')
                }
                tuple(pheno, matched)
            }
    } else {
        // Cartesian product: (pheno chunk) × (genotype chunk), plus cov_pheno.
        step2_input_ch = step1_resolved
            .combine(geno_ch)
            .combine(cov_pheno_ch)
            .map { chunk_id, phenolist, pred, loco, firth_list, firth_gz, covars_used, single_sex, geno, c1, c2, cov_pheno ->
                // Firth gz files may be empty when quant — keep as-is (possibly []).
                def null_files = (firth_gz instanceof List) ? firth_gz : (firth_gz ? [firth_gz] : [])
                tuple(chunk_id, phenolist, pred, loco, firth_list, null_files, cov_pheno, covars_used, single_sex, geno, c1, c2)
            }

        STEP2(step2_input_ch)

        // --- Transpose STEP2 outputs by pheno --------------------------------
        // STEP2 emits a list of <prefix>_<pheno>.regenie.gz files — one per pheno
        // in the chunk. Iterate the phenolist to find each file by suffix (avoids
        // brittle regex on the multi-dot prefix).
        regenie_by_pheno = STEP2.out.main
            .flatMap { chunk_id, phenolist, geno_tag, regenie_files, sex_files, log_files ->
                def regs = (regenie_files instanceof List) ? regenie_files : [regenie_files]
                phenolist.collect { pheno ->
                    def match = regs.find { it.getName().endsWith("_${pheno}.regenie.gz") }
                    tuple(pheno, match)
                }
            }
            .groupTuple()

        sex_by_pheno = STEP2.out.main
            .flatMap { chunk_id, phenolist, geno_tag, regenie_files, sex_files, log_files ->
                def sxs = (sex_files instanceof List) ? sex_files : (sex_files ? [sex_files] : [])
                phenolist.collect { pheno ->
                    def matches = sxs.findAll { f ->
                        def n = f.getName()
                        n.endsWith(".${pheno}.sex_spec.gz") || n.contains('NOT_DONE')
                    }
                    tuple(pheno, matches)
                }
            }
            .groupTuple()
            .map { pheno, lists ->
                def flat = lists.flatten().findAll { it != null }
                tuple(pheno, flat)
            }
    }

    gather_input = regenie_by_pheno
        .join(sex_by_pheno, remainder: true)
        .map { pheno, regs, sexes ->
            tuple(pheno, regs ?: [], sexes ?: [])
        }

    GATHER(gather_input)

    // --- SUMMARY (optional) -----------------------------------------------
    if (!params.skip_summary) {
        fg_annotation_ch = Channel.value(
            tuple(
                file(params.finngen_annotation),
                file("${params.finngen_annotation}.tbi"),
            )
        )

        summary_input = GATHER.out.main
            .map { pheno, reg, sex_diff, sex_diff_tbi, pheweb, pheweb_tbi, qqo, qqe, pngs, qqt ->
                tuple(pheno, pheweb)
            }
            .combine(fg_annotation_ch)
            .map { pheno, pheweb, fg, fg_tbi -> tuple(pheno, pheweb, fg, fg_tbi) }

        SUMMARY(summary_input)

        // --- CODING_GATHER (optional) -------------------------------------
        if (!params.skip_coding_gather) {
            coding_ch = SUMMARY.out.main
                .map { pheno, summary_txt, coding_txt -> coding_txt }
                .collect()
            CODING_GATHER(coding_ch)
        }
    }
}
