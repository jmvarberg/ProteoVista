#visualizations for ProteoVista QC Plots

#Plot 1: Number of Proteins and Peptides per sample
#need UI selection for peptide/protein
#need UI selection for separate or collapsed to groups
ids_by_sample <- function(dataset, type = c("peptide", "protein")) {

    #plot for peptides
    if(type == "peptide") {

        data <- dataset$peptides |>
            dplyr::group_by(sample_id) |>
            dplyr::filter(detect==TRUE) |>
            count() |>
            dplyr::mutate(peptides=n)

        #join in group info
        sample_grouping <- dataset$samples |>
            dplyr::select(sample_id, group) |>
            dplyr::right_join(peptide_data)
    } else if(type == "protein") {

        data <- dataset$proteins |>
            dplyr::group_by(sample_id) |>
            dplyr::filter(detect==TRUE) |>
            count() |>
            dplyr::mutate(peptides=n)

        #join in group info
        sample_grouping <- dataset$samples |>
            dplyr::select(sample_id, group) |>
            dplyr::right_join(peptide_data)

    }

        #if UI input is individual samples
        if(ui$pep_plot_type=="separate") {
            #plot by sample
            plot <- ggplot(sample_grouping, aes(x=peptides, y=sample_id, fill = group)) +
                geom_col() +
                scale_fill_npg(name = "Group") +
                xlab("Number of Detected Peptides") +
                ylab("") +
                theme_bw() +
                theme(legend.position = "top")
        } else if(ui$pep_plot_type=="grouped") {

            #plot by group
            group_summary <- sample_grouping |> dplyr::group_by(group) |> summarise(across(peptides, .fns = list(mean = mean, sd = sd)))

            plot <- ggplot2::ggplot(group_summary, aes(x=round(peptides_mean, 0), y=group, fill = group)) +
                geom_col(width =-0.8) +
                geom_point(data=sample_grouping, aes(x=peptides, y=group), color = "black", size =3, alpha =0.8, show.legend = F) +
                geom_errorbar(data = group_summary, aes(y=group, xmin = peptides_mean - peptides_sd, xmax = peptides_mean + peptides_sd), width = 0.2, linewidth = 1) +
                #geom_text(data = group_summary, aes(y=group, label = round(peptides_mean, 0)), size = 5) +
                geom_label(data = group_summary, aes(label = round(peptides_mean, 0)), fill = "white", color = "black", size =6, nudge_y = 0.25) +
                scale_fill_npg(name = "Group") +
                xlab("Number of Detected Peptides (Group Mean + Std. Dev)") +
                ylab("") +
                theme_cowplot() +
                theme(legend.position = "top")

        }

}


# Value Boxes Functions and UI elements -----------------------------------

output$vbox_total_peptides <- renderText({
    value_box_values()$total_peptides
})

output$vbox_total_proteins <- renderText({
    value_box_values()$total_proteins
})

output$vbox_avg_des <- renderText({
    value_box_values()$avg_de_proteins
})


