#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   dplyr::left_join(readr::read_tsv(file = url_erbb2, col_names = TRUE, col_types = cols(.default = col_character())), by = "sample number") %>%
	   dplyr::select(sample_name = `sample number`,
			 is_tp53_mutated_yes_no = TP53,
			 is_tp53_mutated_0_1 = `TP53 mut?`,
			 hgvsp_short = `IMPACT TP53`,
			 carcinoma_component = `carcinoma component`,
			 carcinoma_classification = `classification of carcinoma (1=serous, 2=endometrioid, 3=HGNOS, 4=undiff)`,
			 sarcoma_classification = `sarcoma component classification (1=homologous; 2=heterologous-including rhabdo, chondroid, chondrosarc)`,
			 `sarcoma_predominant_%` = `sarcoma predominant? >50%`,
			 sarcoma_predominant_yes_no = `sarcoma predominant? (1=yes, 2=no)`,
			 HER2_gene_expression_1 = `HER-2 gene amplification_1`,
			 HER2_gene_expression_2 = `HER-2 gene amplification_2`,
			 HER2_gene_expression_3 = `HER-2 gene amplification_3`,
			 HER2_gene_expression_mean = `HER-2 amplification_avg`,
			 TROP2_gene_expression_1 = `TROP-2 amplification_1`,
			 TROP2_gene_expression_2 = `TROP-2 amplification_2`,
			 TROP2_gene_expression_3 = `TROP-2 amplification_3`,
			 TROP2_gene_expression_mean = `TROP-2 amplification_avg`,
			 ERBB2_gene_amplification) %>%
	   dplyr::mutate(is_tp53_mutated_yes_no = case_when(
		   is_tp53_mutated_yes_no == "Yes" ~ "yes",
		   is_tp53_mutated_yes_no == "No" ~ "no"
	   )) %>%
	   dplyr::mutate(is_tp53_mutated_0_1 = case_when(
		   is_tp53_mutated_0_1 == "1" ~ "1",
		   is_tp53_mutated_0_1 == "2" ~ "0"
	   )) %>%
	   dplyr::mutate(hgvsp_short = gsub(pattern = "(Driver)", replacement = "", x = hgvsp_short, fixed = TRUE)) %>%
	   dplyr::mutate(hgvsp_short = gsub(pattern = ",", replacement = ";", x = hgvsp_short, fixed = TRUE)) %>%
	   dplyr::mutate(hgvsp_short = ifelse(hgvsp_short == "no alteration", NA, hgvsp_short)) %>%
	   dplyr::mutate(carcinoma_classification = case_when(
		   carcinoma_classification == "1" ~ "serous",
		   carcinoma_classification == "2" ~ "endometrioid",
		   carcinoma_classification == "3" ~ "HGNOS",
		   carcinoma_classification == "4" ~ "undifferentiated",
	   )) %>%
	   dplyr::mutate(sarcoma_classification = case_when(
		   is.na(sarcoma_classification) ~ "NA",
		   sarcoma_classification == "1" ~ "homologous",
		   sarcoma_classification == "2" ~ "heterologous"
	   )) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = case_when(
		   grepl("N", `sarcoma_predominant_%`, fixed = TRUE) ~ "0",
		   TRUE ~ `sarcoma_predominant_%`
	   )) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("Y (", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("Y(", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("S:", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("?S:", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("s:", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("C:", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("%", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub(")", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = gsub("?", "", x = `sarcoma_predominant_%`, fixed = TRUE)) %>%
	   dplyr::mutate(`sarcoma_predominant_%` = ifelse(`sarcoma_predominant_%` == "Yes" | `sarcoma_predominant_%` == "Y", NA, `sarcoma_predominant_%`)) %>%
	   dplyr::mutate(sarcoma_predominant_yes_no = case_when(
		   sarcoma_predominant_yes_no == "1" ~ "yes",
		   sarcoma_predominant_yes_no == "2" ~ "no"
	   )) %>%
	   dplyr::mutate(ERBB2_gene_amplification = case_when(
		   is.na(ERBB2_gene_amplification) ~ "Copy-neutral",
		   ERBB2_gene_amplification == "Not Available" ~ "Not available",
		   TRUE ~ ERBB2_gene_amplification
	   )) %>%
	   readr::type_convert()

plot_ = manifest %>%
	reshape2::melt(id.vars = c("sample_name", "ERBB2_gene_amplification"), measure.vars = c("HER2_gene_expression_1", "HER2_gene_expression_2", "HER2_gene_expression_3")) %>%
	dplyr::group_by(sample_name) %>%
	dplyr::summarize(ERBB2_gene_amplification = unique(ERBB2_gene_amplification),
			 mean = mean(value, na.rm = TRUE),
			 min = min(value, na.rm = TRUE),
			 max = max(value, na.rm = TRUE)) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(sample_name = gsub(pattern = "cs", replacement = "CS", x = sample_name, fixed = TRUE)) %>%
	dplyr::mutate(ERBB2_gene_amplification = case_when(
		is.na(ERBB2_gene_amplification) ~ "grey",
		TRUE ~ ERBB2_gene_amplification
	)) %>%
	dplyr::arrange(desc(mean)) %>%
	dplyr::mutate(sample_name = factor(sample_name, levels = unique(sample_name), ordered = TRUE)) %>%
	ggplot(aes(x = sample_name, y = mean, ymin = min, ymax = max, fill = ERBB2_gene_amplification)) +
	geom_bar(stat = "identity", width = .85) +
	geom_pointrange(stat = "identity", shape = 21, fill = "white", color = "grey10", size = .6) +
	scale_fill_brewer(type = "qual", palette = 6) +
	scale_y_sqrt(limits = c(0, .81),
		     breaks = c(seq(from = .01, to = .09, by = .01), seq(from = .1, to = .9, by = .1)),
		     labels = c(c("", ".02", "", ".04", "", ".06", "", ".08", ""), seq(from = .1, to = .9, by = .1))) +
	xlab("") +
	ylab("Expression Fold-Change") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 14),
 	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = guide_legend(title = "ERBB2\nAmplification"))

pdf(file = "../res/Fold_change_by_Sample_Her2.pdf", width = 14, height = 6)
print(plot_)
dev.off()

plot_ = manifest %>%
	reshape2::melt(id.vars = "sample_name", measure.vars = c("TROP2_gene_expression_1", "TROP2_gene_expression_2", "TROP2_gene_expression_3")) %>%
	dplyr::group_by(sample_name) %>%
	dplyr::summarize(mean = mean(value, na.rm = TRUE),
			 min = min(value, na.rm = TRUE),
			 max = max(value, na.rm = TRUE)) %>%
	dplyr::ungroup() %>%
	tidyr::drop_na() %>%
	dplyr::mutate(sample_name = gsub(pattern = "cs", replacement = "CS", x = sample_name, fixed = TRUE)) %>%
	dplyr::arrange(desc(mean)) %>%
	dplyr::mutate(sample_name = factor(sample_name, levels = unique(sample_name), ordered = TRUE)) %>%
	ggplot(aes(x = sample_name, y = mean, ymin = min, ymax = max)) +
	geom_bar(stat = "identity", width = .85) +
	geom_pointrange(stat = "identity", shape = 21, fill = "white", color = "grey10", size = .6) +
	scale_y_sqrt(limits = c(0, 150),
		     breaks = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20)),
		     labels = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20))) +
	xlab("") +
	ylab("Expression Fold-Change") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 14),
 	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Fold_change_by_Sample_Trop2.pdf", width = 14, height = 6)
print(plot_)
dev.off()

