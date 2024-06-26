**Quantifying the global film festival circuit: Networks, diversity, and public value creation**
<br> by Vejune Zemaityte, Andres Karjus, Ulrike Rohn, Maximilian Schich, Indrek Ibrus

Paper here in open access: Zemaityte V, Karjus A, Rohn U, Schich M, Ibrus I (2024) Quantifying the global film festival circuit: Networks, diversity, and public value creation. PLOS ONE 19(3): e0297404. https://doi.org/10.1371/journal.pone.0297404

To replicate the analyses and recreate the graphs, follow these steps:
- Install R (4.2.1 or higher; we recommend using RStudio as the IDE)
- Download the data (csv files) from the persistent repository at https://doi.org/10.6084/m9.figshare.22682794
- This is the sample of data used in the paper, sourced from the Cinando database. We have carefully anonymized all internal Cinando technical ID numbers, and removed identifying information such as names of people and companies.
- Place all files from figshare in a folder of choice along with the R script files (or the folder reported by `getwd()` in R). 
- Open the `quantifying_festivalcircuit.R` file in R or RStudio
- Follow the steps in `quantifying_festivalcircuit.R`: source the functions R script file (this also installs packages if missing), and run all lines in the file to recreate the analyses and graphs; the graphs will be exported into the working directory.

If you make use of these resources, please cite the paper:<br>
Vejune Zemaityte, Andres Karjus, Ulrike Rohn, Maximilian Schich, Indrek Ibrus (2023). "Quantifying the global film festival circuit: Networks, diversity, and public value creation". Socarxiv preprint.
<br>(the paper is currently under review; reference will be updated once published)

A dashboard with interactive versions of the graphs from the paper is available at: https://andreskarjus.github.io/cinandofestivals
