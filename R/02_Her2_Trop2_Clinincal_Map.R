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
	   dplyr::left_join(readr::read_tsv(file = url_ihc, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			    readr::type_convert() %>%
			    dplyr::rename(sample_name = `sample number`,
					  her2_ihc = `HER-2 IHC score`,
					  trop2_ihc = `TROP-2 IHC score`), by = "sample_name")

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
	dplyr::ungroup() %>%
	dplyr::mutate(sample_name = gsub(pattern = "cs", replacement = "CS", x = sample_name, fixed = TRUE)) %>%
	dplyr::arrange(desc(her2))

pdf(file = "../res/Her2_Strip_Chart.pdf", width = 14, height = 6)
draw(Heatmap(matrix = smry_ %>%
	     	      dplyr::select(-her2) %>%
	     	      as.data.frame() %>%
	     	      tibble::column_to_rownames(var = "sample_name") %>%
	     	      t(),
	     col = c("1" = "#fb9a99", "2" = "#1f78b4", "3" = "#b2df8a", "4" = "#33a02c", "100" = "#bae4bc", "101" = "#7bccc4", "102" = "#43a2ca", "103" = "#0868ac"),
	     name = " ",
	     na_col = "#f0f0f0",
	     rect_gp = gpar(col = "white"),
	     border = NA,
	     
	     cluster_rows = FALSE,
	     show_row_dend = FALSE,
	     
	     row_names_side = "left",
	     row_names_gp = gpar(fontsize = 9),
	     show_row_names = TRUE,
	     
	     column_names_side = "bottom",
	     column_names_gp = gpar(fontsize = 9),
	     show_column_names = TRUE,
	     
	     cluster_columns = FALSE,
	     column_order = smry_ %>% .[["sample_name"]],
	     use_raster = FALSE,
	     width = unit(14*2, "cm"),
	     height = unit(2, "cm"),

	     show_heatmap_legend = TRUE,
	     heatmap_legend_param = list(legend_height = unit(3, "cm"), legend_width = unit(2, "cm"))))
dev.off()

smry_ = manifest %>%
	dplyr::filter(!(sample_name %in% exclude)) %>%
	reshape2::melt(id.vars = c("sample_name", "tp53_mutated", "carcinoma_classification", "sarcoma_classification", "sarcoma_predominant", "trop2_ihc"),
		       measure.vars = c("TROP2_gene_expression_1", "TROP2_gene_expression_2", "TROP2_gene_expression_3")) %>%
	dplyr::group_by(sample_name) %>%
	dplyr::summarize(trop2 = mean(value, na.rm = TRUE),
			 tp53_mutated = unique(tp53_mutated),
			 carcinoma_classification = unique(carcinoma_classification),
			 sarcoma_classification = unique(sarcoma_classification),
			 sarcoma_predominant = unique(sarcoma_predominant),
			 trop2_ihc = mean(trop2_ihc, na.rm=TRUE)) %>%
	dplyr::ungroup() %>%
	tidyr::drop_na() %>%
	dplyr::mutate(sample_name = gsub(pattern = "cs", replacement = "CS", x = sample_name, fixed = TRUE)) %>%
	dplyr::arrange(trop2_ihc) %>%
	dplyr::mutate(trop2_cat = case_when(
		trop2_ihc >= 0 & trop2_ihc < 50 ~ "100",
		trop2_ihc >= 50 & trop2_ihc < 100 ~ "101",
		trop2_ihc >= 100 & trop2_ihc < 150 ~ "102",
		trop2_ihc >= 150 & trop2_ihc < 200 ~ "103",
		trop2_ihc >= 200 & trop2_ihc < 250 ~ "104",
		trop2_ihc >= 250 & trop2_ihc <= 300 ~ "105",
		TRUE ~ "NA"
	))

palette = colorRampPalette(colors = c("#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"))(length(unique(smry_ %>% .[["trop2_cat"]])))
names(palette) = sort(unique(smry_ %>% .[["trop2_cat"]]))

pdf(file = "../res/Trop2_Strip_Chart.pdf", width = 14, height = 6)
draw(Heatmap(matrix = smry_ %>%
	     	      dplyr::select(-trop2) %>%
	     	      dplyr::select(-trop2_ihc) %>%
	     	      as.data.frame() %>%
	     	      tibble::column_to_rownames(var = "sample_name") %>%
	     	      t(),
	     col = c("1" = "#fb9a99", "2" = "#1f78b4", "3" = "#b2df8a", "4" = "#33a02c", "NA" = "dedbd8", palette),
	     name = " ",
	     na_col = "#f0f0f0",
	     rect_gp = gpar(col = "white"),
	     border = NA,
	     
	     cluster_rows = FALSE,
	     show_row_dend = FALSE,
	     
	     row_names_side = "left",
	     row_names_gp = gpar(fontsize = 9),
	     show_row_names = TRUE,
	     
	     column_names_side = "bottom",
	     column_names_gp = gpar(fontsize = 9),
	     show_column_names = TRUE,
	     
	     cluster_columns = FALSE,
	     column_order = smry_ %>% .[["sample_name"]],
	     use_raster = FALSE,
	     width = unit(14*2, "cm"),
	     height = unit(2, "cm"),

	     show_heatmap_legend = TRUE,
	     heatmap_legend_param = list(legend_height = unit(3, "cm"), legend_width = unit(2, "cm"))))
dev.off()
