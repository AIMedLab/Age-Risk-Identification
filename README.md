# Age-Risk-Identification
## 1. Introduction
This repository contains source code for paper "A Computational Framework for Identifying Age Risks in Drug-Adverse Event Pairs".

In this paper, we developed a statistical computational framework to detect the age group of patients who are susceptible to some ADEs after taking specific drugs. We applied our matheodology to FDA Adverse Event Reporting System (FAERS) and indentified age-associated drug-adverse event pairs as well as their highest age risk group.

We divided patients into four age groups, 0-14 years old as children, 15-24 years old as youth, 25-64 years old as adult and \textgreater 65 years old as senior on the basis of the criteria provided by World Population Prospects from The United Nations Department of Economic and Social Affairs(available at https://population.un.org/wpp/).

The age risks are detected mainly through four steps: discover age differences for specific drugs, discover age differences for specific drug-ADE pairs, discover the age group with higher risk, remove confounding effect by age bias. The illustration of our framework can be found in Fig.1.

<img width="671" alt="Screen Shot 2021-03-01 at 9 44 48 PM" src="https://user-images.githubusercontent.com/79823323/110710818-8800d700-81cc-11eb-888c-2a79d36e6643.png">

Fig.1:Illustration of the proposed methodology: In the first step, overall Chi-squared tests for each drug are performed to identify drugs with overall age difference. Then overall Chi-squared tests for each identified drugs and ADEs pair are performed to identify drug-ADE pairs with age difference. Next Chi-squared tests for age group comparisons within each pair are performed. RORs for every two age groups are computed and ranked, which quantifies the age risks. At the end, a logistic model is built for each detected pair and the Likelihood Ratio Test is performed on the interaction of drug and age group to remove age bias.

## 2. Dataset
FAERS: a database that contains information on adverse event and medication error reports submitted to FDA. We collected FAERS quarterly submissions from 2004 to the third quarter of 2018 and cleaned and normalized the data(Banda, Juan M. et al., 2017).



