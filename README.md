**Quantifying the global film festival circuit: Networks, diversity, and public value creation**
<br> by Vejune Zemaityte, Andres Karjus, Ulrike Rohn, Maximilian Schich, Indrek Ibrus

Preprint available here: [link]

To replicate the analyses and recreate the graphs, follow these steps:
- Install R (4.2.1 or higher; we recommend using RStudio as the IDE)
- Unpack the data (csv files) from the zip file and place them in a folder of choice along with the script files (or the folder reported by `getwd()` in R). This is the sample of data used in the paper, sourced from the Cinando database. We have carefully anonymized all internal Cinando technical ID numbers, and removed identifying information such as names of people and companies.
- Open the `quantifying_festivalcircuit.R` file in R or RStudio
- Follow the steps in `quantifying_festivalcircuit.R`: source the functions R script file (this also installs packages if missing), and run all lines in the file to recreate the analyses and graphs; the graphs will be exported into the working directory.

If you make use of these resources, please cite the paper above.

A dashboard with interactive versions of the graphs of the paper is available here: https://andreskarjus.github.io/cinandofestivals
