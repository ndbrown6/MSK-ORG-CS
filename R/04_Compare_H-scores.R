#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

smry_ = readr::read_tsv(file = url_trop2, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	dplyr::filter(`Carcinoma H-score` != "n/a") %>%
	dplyr::filter(`Overall H-score` != "n/a") %>%
	dplyr::filter(`% Sarcoma` != "n/a") %>%
	readr::type_convert() %>%
	tidyr::drop_na() %>%
	dplyr::mutate(`Sarcoma H-score` = 0)


plot_ = smry_ %>%
	dplyr::select(`sample number`, `Carcinoma H-score`, `Sarcoma H-score`) %>%
	reshape2::melt() %>%
	dplyr::mutate(variable = gsub(" H-score", "", variable, fixed = TRUE)) %>%
	ggplot(aes(x = variable, y = value)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	scale_x_discrete() +
	scale_y_continuous() +
	xlab(" ") +
	ylab("TROP2 H-score") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = 300,
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("Carcinoma", "Sarcoma"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11))
	
pdf(file = "../res/Trop2_H-score.pdf", width = 3, height = 4)
print(plot_)
dev.off()

plot_ = smry_ %>%
	dplyr::select(`sample number`, `Carcinoma H-score`, `Sarcoma H-score`) %>%
	reshape2::melt() %>%
	dplyr::mutate(variable = gsub(" H-score", "", variable, fixed = TRUE)) %>%
	ggplot(aes(x = variable, y = value)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_path(data = smry_ %>%
		  	 dplyr::select(`sample number`, `Carcinoma H-score`, `Sarcoma H-score`) %>%
		  	 reshape2::melt() %>%
		  	 dplyr::mutate(variable = gsub(" H-score", "", variable, fixed = TRUE)),
		  mapping = aes(x = variable, y = value, group = `sample number`),
		  size = .1) +
	scale_x_discrete() +
	scale_y_continuous() +
	xlab(" ") +
	ylab("TROP2 H-score") +
	geom_signif(stat = "signif", test = "wilcox.test", test.args = list(exact = FALSE, alternative = "two.sided"),
		    y_position = 300,
		    textsize = 3,
		    vjust = 0,
		    comparison = list(c("Carcinoma", "Sarcoma"))) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 11),
	      axis.text.y = element_text(size = 11))
	
pdf(file = "../res/Trop2_H-score_with_lines.pdf", width = 3, height = 4)
print(plot_)
dev.off()
