# population-genetics-evolution
Simulation of the evolution of a population of N diploid individuals (e.g  humans), study how allele frequencies change over time and how this is impacted by selection

A program that simulates the evolution of N diploid individuals with L SNPs.\
In this simulation, each individual is represented as a pair of 
binary vectors (0=reference allele, 1=alternate allele), corresponding to the alleles on 
the maternal and paternal chromosome.\
For example, if L=5, an individual could be represented as:\
Maternal chrom: 0 0 1 0 1\
Paternal chrom: 0 1 1 0 0\
This would mean that this individual is homozygous for the reference allele at 
positions 0 and 3, heterozygous at positions 1 and 4, and homozygous for the 
alternate allele at position 2.\
\
The simulation starts from some initial population (founding members), and then 
generates a next generation by repeatedly choosing two individuals as parents and 
producing an offspring from them. We will work with the following assumptions:
- Initial population:
    - The phased genotype information for the initial population is stored in 
the file initial_population.vcf. Note for each SNP, exactly one individual is 
heterozygous, while the others are homozygous for the reference allele. 
This means that the alternate allele frequencies of every SNP is 1/(2N). This mimics an initial population with no standing variation, except for 100 mutations specific to each individual.
- Evolution:
    - We will simulate a panmixing population (i.e. a population where any two 
individuals can mate with any other; there are no males or females) of 
constant size N.
    - To go from a population at generation G to a population at generation 
G+1:\
        &nbsp;&nbsp;&nbsp;&nbsp;For i in 0…N-1\
          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pick parent P1 based on fitness wrt population Pop<sub>G</sub>\
          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pick parent P2≠P1 based on fitness wrt population Pop<sub>G</sub>\
          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Child<sub>i</sub> = reproduce(P1,P2)\
       &nbsp;&nbsp;&nbsp;&nbsp;Pop<sub>G+1</sub> = {Child<sub>0</sub>, …, Child<sub>N-1</sub>}
    - Parent Selection:
        - At each of the N iterations of the loop above, the probability that 
individual i is selected as parent P<sub>1</sub> is  $$fitness(i)/\sum_{j=1}^N fitness(j)$$
        - The probability that individual i′ ≠ P<sub>1</sub> is selected as P<sub>2</sub> is  $$fitness(i)/\sum_{j=1...N;j≠P1} fitness(j)$$
        - Note: A given individual can be selected as parent any number of times. However, an individual cannot at the same time be the mother and the father
    - Reproduction:\
       Suppose we have:\
       &nbsp;&nbsp;&nbsp;&nbsp;Parent P<sub>1</sub> (mother) with maternal chrom a<sub>0</sub>,…,a<sub>L-1</sub> and paternal chrom b<sub>0</sub>,…,b<sub>L-1</sub>\
      &nbsp;&nbsp;&nbsp;&nbsp;Parent P<sub>2</sub> (father) with maternal chrom c<sub>0</sub>,…,c<sub>L-1</sub> and paternal chrom d<sub>0</sub>,…,d<sub>L-1</sub>
      - Choose crossing-over position k uniformly among {0, 1, …, L-1}
      - Child maternal chromosome = a<sub>0</sub>…a<sub>k</sub>b<sub>k+1</sub>…b<sub>L-1</sub> or b<sub>0</sub>…b<sub>k</sub>a<sub>k+1</sub>…a<sub>L-1</sub> (with probability 0.5 each)
      - Choose crossing-over posi9on k’ uniformly among {0, 1, …, L-1}
      -  Child paternal chromosome = c<sub>0</sub>…c<sub>k'</sub>d<sub>k'+1</sub>…d<sub>L-1</sub> or d<sub>0</sub>…d<sub>k'</sub>c<sub>k'+1</sub>…c<sub>L-1</sub> (with probability 0.5 each)
    - Length of the simulation: Code is be able to simulate the 
evolution of an initial population over any user-provided number of 
generations.

### Simulation of the evolutionof over 20 generations, assuming alleles are selectively neutral
![image](https://github.com/tgompper1/population-genetics-evolution/assets/33359970/ba873ecc-b638-49d7-8331-63fcf01e6f92)

