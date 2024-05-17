# TLC to Chromatogram simulator

![](https://github.com/gilbertblanson/TLC-to-Chromatogram/blob/main/demoapp.gif)

An interactive ShinyApp that takes your TLC results and turns them into predicted chromatograms for manual silica column chromatography.

## Getting Started

You can access the app on any device at [ansonc.shinyapps.io/TLC-Chromatogram/](ansonc.shinyapps.io/TLC-Chromatogram/). To run it locally, download the app.R and run it in RStudio.

Words of warning when using the program:
-  This was specifically developed for manual flash columns using normal phase silica, 40–63 µm (230–400 mesh), 60 Å pore size, pH range of 6.5–7.5. The transference from a TLC Rf to a normal phase elution volume (V'r) uses an empirical correction constant that was determined by Justin Fair and Chad Kormos ([link](https://doi.org/10.1016/j.chroma.2008.09.085)). Commercial silica cartridges have different efficiencies and thus the observed V'r and bandwidth may differ.
-  As per the original publication: Greater error for more highly retained analytes stems from the inverse relationship between Rf and V'r and thus, errors are more likely to propagate. It is therefore recommended that the spreadsheet not be used for Rf values less than ∼0.08; RF values less than ∼0.1 produce extremely broad peak shapes due to the broadening associated with a large V'r. Due to the nature in which bandwidth (Vb) and resolution (Rs) are calculated, any error in the calculation of V'r will be carried into both Vb and Rs.

## Acknowledgments

-   Justin Fair and Chad Kormos for their [original publication](https://doi.org/10.1016/j.chroma.2008.09.085) which presented the theory this ShinyApp builds upon.
-   Pavel Jandera and Jaroslav Churáček for their work on [capacity factor prediction](https://doi.org/10.1016/S0021-9673(00)99325-7), and Paweł Kręcisz, Kamila Czarnecka, and Paweł Szymański for [their work](https://doi.org/10.1093/chromsci/bmab097) which provided the inspiration to translate that to Rf prediction.
-   MKB for his assistance and guidance throughout the project.
