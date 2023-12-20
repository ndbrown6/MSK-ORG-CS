#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   dplyr::select(sample_name = `sample number`,
			 is_tp53_mutated_yes_no = TP53,
			 is_tp53_mutated_0_1 = `TP53 mut?`,
			 hgvsp_short = `IMPACT TP53`,
			 carcinoma_component = `carcinoma component`,
			 carcinoma_classification = `classification of carcinoma (1=serous, 2=endometrioid, 3=HGNOS, 4=undiff)`,
			 sarcoma_classification = `sarcoma component classification (1=homologous; 2=heterologous-including rhabdo, chondroid, chondrosarc)`,
			 `sarcoma_predominant_%` = `sarcoma predominant? >50%`,
			 sarcoma_predominant_yes_no = `sarcoma predominant? (1=yes, 2=no)`,
			 her2_amplification_rep_1 = `HER-2 gene amplification`,
			 her2_amplification_rep_2 = X12,
			 her2_amplification_rep_3 = X13,
			 her2_amplification_average = `HER-2 amplification avg`,
			 trop2_amplification_rep_1 = `TROP-2 amplification`,
			 trop2_amplification_rep_2 = X16,
			 trop2_amplification_rep_3 = X17,
			 trop2_amplification_average = `TROP-2 amplification average`) %>%
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
		   carcinoma_classification == "1" ~ "homologous",
		   carcinoma_classification == "2" ~ "heterologous"
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
	   dplyr::mutate(trop2_amplification_average = ifelse(trop2_amplification_average == "#DIV/0!", NA, trop2_amplification_average)) %>%
	   readr::type_convert()
	   

plot_ = manifest %>%
	dplyr::group_by(sample_name) %>%
	dplyr::summarize(her2_amplification_average = her2_amplification_average,
			 her2_amplification_min = min(her2_amplification_rep_1, her2_amplification_rep_2, her2_amplification_rep_3, na.rm=TRUE),
			 her2_amplification_max = max(her2_amplification_rep_1, her2_amplification_rep_2, her2_amplification_rep_3, na.rm=TRUE)) %>%
	dplyr::arrange(desc(her2_amplification_average)) %>%
	dplyr::mutate(sample_name = factor(sample_name, levels = unique(sample_name), ordered = TRUE)) %>%
	ggplot(aes(x = sample_name, y = her2_amplification_average, ymin = her2_amplification_min, ymax = her2_amplification_max)) +
	geom_bar(stat = "identity", width = .85) +
	geom_pointrange(stat = "identity", shape = 21, fill = "white", color = "grey10", size = .65) +
	scale_y_sqrt(limits = c(0, .81),
		     breaks = c(seq(from = .01, to = .09, by = .01), seq(from = .1, to = .9, by = .1)),
		     labels = c(c("", ".02", "", ".04", "", ".06", "", ".08", ""), seq(from = .1, to = .9, by = .1))) +
	xlab("") +
	ylab("Expression Fold-Change") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 14),
 	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Fold_change_by_Sample.pdf", width = 14, height = 6)
print(plot_)
dev.off()
