#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   dplyr::select(sample_name = `sample number`,
			 tp53_mutated = `TP53 mut?`,
			 hgvsp_short = `IMPACT TP53`,
			 carcinoma_component = `carcinoma component`,
			 carcinoma_classification = `classification of carcinoma (1=serous, 2=endometrioid, 3=HGNOS, 4=undiff)`,
			 sarcoma_classification = `sarcoma component classification (1=homologous; 2=heterologous-including rhabdo, chondroid, chondrosarc)`,
			 `sarcoma_predominant_%` = `sarcoma predominant? >50%`,
			 sarcoma_predominant = `sarcoma predominant? (1=yes, 2=no)`,
			 HER2_gene_expression_1,
			 HER2_gene_expression_2,
			 HER2_gene_expression_3,
			 HER2_gene_expression_mean,
			 TROP2_gene_expression_1,
			 TROP2_gene_expression_2,
			 TROP2_gene_expression_3,
			 TROP2_gene_expression_mean,
			 ERBB2_gene_amplification) %>%
	   dplyr::mutate(hgvsp_short = gsub(pattern = "(Driver)", replacement = "", x = hgvsp_short, fixed = TRUE)) %>%
	   dplyr::mutate(hgvsp_short = gsub(pattern = ",", replacement = ";", x = hgvsp_short, fixed = TRUE)) %>%
	   dplyr::mutate(hgvsp_short = ifelse(hgvsp_short == "no alteration", NA, hgvsp_short)) %>%
	   readr::type_convert() %>%
	   dplyr::left_join(readr::read_tsv(file = url_ihc_her2, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			    readr::type_convert() %>%
			    dplyr::rename(sample_name = `sample number`,
					  her2_ihc = `HER-2 IHC score`),
			    by = "sample_name") %>%
  	   dplyr::left_join(readr::read_tsv(file = url_ihc_trop2, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			    readr::type_convert() %>%
			    dplyr::rename(sample_name = `sample number`,
					  trop2_carcinoma_ihc = `TROP2 carcinoma H score`,
					  trop2_overall_ihc = `TROP 2 overall H score`),
			    by = "sample_name")

smry_ = manifest %>%
	reshape2::melt(id.vars = c("sample_name", "tp53_mutated", "carcinoma_classification", "sarcoma_classification", "sarcoma_predominant", "her2_ihc"),
		       measure.vars = c("HER2_gene_expression_1", "HER2_gene_expression_2", "HER2_gene_expression_3")) %>%
	dplyr::group_by(sample_name) %>%
	dplyr::summarize(her2 = mean(value, na.rm = TRUE),
			 tp53_mutated = unique(tp53_mutated),
			 carcinoma_classification = unique(carcinoma_classification),
			 sarcoma_classification = unique(sarcoma_classification),
			 sarcoma_predominant = unique(sarcoma_predominant),
			 her2_ihc = unique(her2_ihc)+100) %>%
	dplyr::ungroup()
	
plot_ = smry_ %>%
        dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff",
	)) %>%
        dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_predominant, y = her2)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_classification, y = her2, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Yes", "No")) +
	scale_y_sqrt(limits = c(0, .81),
		     breaks = c(seq(from = .01, to = .09, by = .01), seq(from = .1, to = .9, by = .1)),
		     labels = c(c("", ".02", "", ".04", "", ".06", "", ".08", ""), seq(from = .1, to = .9, by = .1))) +
	xlab("Predominant sarcoma") +
	ylab("HER2 Expression Fold-Change") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = sqrt(0.75),
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Her2_by_Predominat_Sarcoma.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff",
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_predominant) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Her2_by_Predominat_Sarcoma")

plot_ = smry_ %>%
	dplyr::filter(!is.na(sarcoma_classification)) %>%
        dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_classification, y = her2)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_classification, y = her2, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Homologous", "Heterologous")) +
	scale_y_sqrt(limits = c(0, .81),
		     breaks = c(seq(from = .01, to = .09, by = .01), seq(from = .1, to = .9, by = .1)),
		     labels = c(c("", ".02", "", ".04", "", ".06", "", ".08", ""), seq(from = .1, to = .9, by = .1))) +
	xlab("Sarcoma classification") +
	ylab("HER2 Expression Fold-Change") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = sqrt(0.75),
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Her2_by_Sarcoma_Classification.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::filter(!is.na(sarcoma_classification)) %>%
dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_classification) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Her2_by_Sarcoma_Classification")


smry_ = manifest %>%
	dplyr::filter(!(sample_name %in% exclude)) %>%
	reshape2::melt(id.vars = c("sample_name", "tp53_mutated", "carcinoma_classification", "sarcoma_classification", "sarcoma_predominant", "trop2_overall_ihc", "trop2_carcinoma_ihc"),
		       measure.vars = c("TROP2_gene_expression_1", "TROP2_gene_expression_2", "TROP2_gene_expression_3")) %>%
	dplyr::group_by(sample_name) %>%
	dplyr::summarize(trop2 = mean(value, na.rm = TRUE),
			 tp53_mutated = unique(tp53_mutated),
			 carcinoma_classification = unique(carcinoma_classification),
			 sarcoma_classification = unique(sarcoma_classification),
			 sarcoma_predominant = unique(sarcoma_predominant),
			 trop2_overall_ihc = unique(trop2_overall_ihc),
			 trop2_carcinoma_ihc = unique(trop2_carcinoma_ihc)) %>%
	dplyr::ungroup()
	
plot_ = smry_ %>%
        dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_predominant, y = trop2)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_predominant, y = trop2, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Yes", "No")) +
	scale_y_sqrt(limits = c(0, 150),
		     breaks = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20)),
		     labels = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20))) +
	xlab("Predominant sarcoma") +
	ylab("TROP2 Expression Fold-Change") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = sqrt(140),
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_by_Predominat_Sarcoma.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_predominant) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_by_Predominat_Sarcoma")


plot_ = smry_ %>%
	dplyr::filter(!is.na(sarcoma_classification)) %>%
	dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
		dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_classification, y = trop2)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_classification, y = trop2, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Homologous", "Heterologous")) +
	scale_y_sqrt(limits = c(0, 150),
		     breaks = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20)),
		     labels = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20))) +
	xlab("Sarcoma classification") +
	ylab("TROP2 Expression Fold-Change") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = sqrt(140),
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_by_Sarcoma_Classification.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::filter(!is.na(sarcoma_classification)) %>%
dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_classification) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_by_Sarcoma_Classification")


plot_ = smry_ %>%
	dplyr::filter(!is.na(sarcoma_classification)) %>%
	dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
		dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_predominant:sarcoma_classification, y = trop2)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_predominant:sarcoma_classification, y = trop2, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1:1", "1:2", "2:1", "2:2"),
			 labels = c("Pred sarc\nHom sarc", "Pred sarc\nHet sarc",
				    "Pred car\nHom sarc", "Pred car\nHet sarc")) +
	scale_y_sqrt(limits = c(0, 300),
		     breaks = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20)),
		     labels = c(seq(from = 1, to = 7, by = 2), seq(from = 10, to = 150, by = 20))) +
	xlab("Sarcoma classification") +
	ylab("TROP2 Expression Fold-Change") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = sqrt(seq(140, 300, length = 6)),
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1:1", "1:2"),
				      c("1:1", "2:1"),
				      c("1:1", "2:2"),
				      c("1:2", "2:1"),
				      c("1:2", "2:2"),
				      c("2:1", "2:2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_by_Predominant_Sarcome,Sarcoma_Classification.pdf", width = 5.5, height = 5.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::filter(!is.na(sarcoma_classification)) %>%
dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_predominant:sarcoma_classification) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_by_Predominant_Sarcome,Sarcoma_Classification")


plot_ = smry_ %>%
        dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_predominant, y = trop2_overall_ihc)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_predominant, y = trop2_overall_ihc, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Yes", "No")) +
	scale_y_continuous(limits = c(0, 350),
			   breaks = c(0, 50, 100, 150, 200, 250, 300),
			   labels = c(0, 50, 100, 150, 200, 250, 300)) +
	xlab("Predominant sarcoma") +
	ylab("TROP2 H-score") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = 300,
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_IHC_by_Predominat_Sarcoma.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_predominant) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_IHC_by_Predominat_Sarcoma")


plot_ = smry_ %>%
	dplyr::filter(!is.na(sarcoma_classification)) %>%
	dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
		dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_classification, y = trop2_overall_ihc)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_classification, y = trop2_overall_ihc, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Homologous", "Heterologous")) +
	scale_y_continuous(limits = c(0, 350),
			   breaks = c(0, 50, 100, 150, 200, 250, 300),
			   labels = c(0, 50, 100, 150, 200, 250, 300)) +
	xlab("Sarcoma classification") +
	ylab("TROP2 H-score") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = 300,
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_IHC_by_Sarcoma_Classification.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::filter(!is.na(sarcoma_classification)) %>%
dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_classification) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_IHC_by_Sarcoma_Classification")


plot_ = smry_ %>%
	dplyr::filter(!is.na(sarcoma_classification)) %>%
	dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
		dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_predominant:sarcoma_classification, y = trop2_overall_ihc)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_predominant:sarcoma_classification, y = trop2_overall_ihc, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1:1", "1:2", "2:1", "2:2"),
			 labels = c("Pred sarc\nHom sarc", "Pred sarc\nHet sarc",
				    "Pred car\nHom sarc", "Pred car\nHet sarc")) +
	scale_y_continuous(limits = c(0, 350),
			   breaks = c(0, 50, 100, 150, 200, 250, 300),
			   labels = c(0, 50, 100, 150, 200, 250, 300)) +
	xlab("Sarcoma classification") +
	ylab("TROP2 H-score") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = seq(265, 350, length = 6),
		    textsize = 3,
		    vjust = 0,
		    tip_length = 0.01,
		    comparison = list(c("1:1", "1:2"),
				      c("1:1", "2:1"),
				      c("1:1", "2:2"),
				      c("1:2", "2:1"),
				      c("1:2", "2:2"),
				      c("2:1", "2:2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_IHC_by_Predominant_Sarcome,Sarcoma_Classification.pdf", width = 5.5, height = 5.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::filter(!is.na(sarcoma_classification)) %>%
dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_predominant:sarcoma_classification) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_IHC_by_Predominant_Sarcome,Sarcoma_Classification")


plot_ = smry_ %>%
        dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_predominant, y = trop2_carcinoma_ihc)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_predominant, y = trop2_carcinoma_ihc, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Yes", "No")) +
	scale_y_continuous(limits = c(0, 350),
			   breaks = c(0, 50, 100, 150, 200, 250, 300),
			   labels = c(0, 50, 100, 150, 200, 250, 300)) +
	xlab("Predominant sarcoma") +
	ylab("TROP2 H-score (carcinoma)") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = 300,
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_IHC_Carcinoma_by_Predominat_Sarcoma.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_predominant) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_IHC_Carcinoma_by_Predominat_Sarcoma")

plot_ = smry_ %>%
	dplyr::filter(!is.na(sarcoma_classification)) %>%
	dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
		dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_classification, y = trop2_carcinoma_ihc)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_classification, y = trop2_carcinoma_ihc, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1", "2"),
			 labels = c("Homologous", "Heterologous")) +
	scale_y_continuous(limits = c(0, 350),
			   breaks = c(0, 50, 100, 150, 200, 250, 300),
			   labels = c(0, 50, 100, 150, 200, 250, 300)) +
	xlab("Sarcoma classification") +
	ylab("TROP2 H-score (carcinoma)") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = 300,
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("1", "2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_IHC_Carcinoma_by_Sarcoma_Classification.pdf", width = 4.5, height = 4.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::filter(!is.na(sarcoma_classification)) %>%
dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_classification) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_IHC_Carcinoma_by_Sarcoma_Classification")


plot_ = smry_ %>%
	dplyr::filter(!is.na(sarcoma_classification)) %>%
	dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
	dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
		dplyr::mutate(carcinoma_classification = case_when(
		carcinoma_classification == 1 ~ "Serous",
		carcinoma_classification == 2 ~ "Endometrioid",
		carcinoma_classification == 3 ~ "HGNOS",
		carcinoma_classification == 4 ~ "Undiff"
	)) %>%
	dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
	ggplot(aes(x = sarcoma_predominant:sarcoma_classification, y = trop2_carcinoma_ihc)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(mapping = aes(x = sarcoma_predominant:sarcoma_classification, y = trop2_carcinoma_ihc, color = carcinoma_classification),
		    position = position_jitter(0.15, 0), size = 3.5, shape = 21, fill = "white", alpha = .75, inherit.aes = TRUE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_discrete(breaks = c("1:1", "1:2", "2:1", "2:2"),
			 labels = c("Pred sarc\nHom sarc", "Pred sarc\nHet sarc",
				    "Pred car\nHom sarc", "Pred car\nHet sarc")) +
	scale_y_continuous(limits = c(0, 350),
			   breaks = c(0, 50, 100, 150, 200, 250, 300),
			   labels = c(0, 50, 100, 150, 200, 250, 300)) +
	xlab("Sarcoma classification") +
	ylab("TROP2 H-score (carcinoma)") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = seq(265, 350, length = 6),
		    textsize = 3,
		    vjust = 0,
		    tip_length = 0.01,
		    comparison = list(c("1:1", "1:2"),
				      c("1:1", "2:1"),
				      c("1:1", "2:2"),
				      c("1:2", "2:1"),
				      c("1:2", "2:2"),
				      c("2:1", "2:2"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11)) +
	guides(color = guide_legend(title = "Carcinoma\nclassification"))
	
pdf(file = "../res/Trop2_IHC_Carcinoma_by_Predominant_Sarcome,Sarcoma_Classification.pdf", width = 5.5, height = 5.5)
print(plot_)
dev.off()

smry_ %>%
dplyr::filter(!is.na(sarcoma_classification)) %>%
dplyr::mutate(sarcoma_classification = factor(sarcoma_classification, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(sarcoma_predominant = factor(sarcoma_predominant, levels = c(1,2), ordered = TRUE)) %>%
dplyr::mutate(carcinoma_classification = case_when(
	carcinoma_classification == 1 ~ "Serous",
	carcinoma_classification == 2 ~ "Endometrioid",
	carcinoma_classification == 3 ~ "HGNOS",
	carcinoma_classification == 4 ~ "Undiff"
)) %>%
dplyr::mutate(carcinoma_classification = factor(carcinoma_classification, levels = c("Serous", "Endometrioid", "HGNOS", "Undiff", "NA"), ordered = TRUE)) %>%
dplyr::group_by(sarcoma_predominant:sarcoma_classification) %>%
dplyr::summarize(n = n()) %>%
pander::pander(caption = "Trop2_IHC_Carcinoma_by_Predominant_Sarcome,Sarcoma_Classification")