# Tess Gompper #260947251
import pandas as pd
from random import randrange, randint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
from array import *

FILENAME = './initial_population.vcf'
NUM_SNPS = 10000
NUM_INDIVIDUALS = 100

def population_evolution(N, L, g, fitness_func):
    """ 
        Simulate Population Evolution.
        N = population size 
        L = number of snps
        g = number of generations to be evolved
        fitness_func = provided fitness function, 
            takes population as input, reurns an int representing most fit 
            individual
    """
    # Initial Population
    #   This mimics an initial population with no standing variation, except for
    #   100 mutations specific to each individual.
    initial_pop = pd.read_csv(FILENAME, sep='\t')
    initial_pop = initial_pop.drop(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL'], axis=1)
   
    # Evolution
    #   We will simulate a panmixing population (i.e. a population where any two
    #   individuals can mate with any other; there are no males or females) of 
    #   constant size N.
    parent_population = initial_pop
    # go from a population at generation G to a population at generation G+1:
    for j in range(0, g):
        child_population = evolve(N, L, parent_population, fitness_func)
        parent_population = child_population
    return child_population

def evolve(N, L, parent_population, fitness_func):
    # go from population G to G+1
    child_population = pd.DataFrame()
    for i in range(0, N):
        # pick parent P1 based on fitness wrt population PopG
        p1 = fitness_func(parent_population)

        # pick parent P2 != P1 based on fitness wrt population PopG (done in fitness_func)
        p2 = fitness_func(parent_population, p1)

        str1 = "IND"+str(p1)
        str2 = "IND"+str(p2)
        mother = parent_population[str1].to_list()
        father = parent_population[str2].to_list()

        # Reproduction
        child, child_name = reproduce(mother, father, L, i)
        child_population[child_name] = child

    return child_population


def reproduce(mother, father, L, num):
    # Parent P1 (mother) with maternal chrom a0,…,aL-1 and paternal chrom b0,…,bL-1
    p1=[]
    for snp in mother:
        s = snp.split('|') 
        s[0] = int(s[0])
        s[1] = int(s[1])
        p1.append(s)
        
    # Parent P2 (father) with maternal chrom c0,…,cL-1 and paternal chrom d0,…,dL-1
    p2=[]
    for snp in father:
        s = snp.split('|') 
        s[0] = int(s[0])
        s[1] = int(s[1])
        p2.append(s)
    # Choose crossing-over position k uniformly among {0, 1, …, L-1}
    k = randint(0, L-1)
    # Child maternal chromosome = a0…ak bk+1…bL-1 or b0…bk ak+1…aL-1 (with probability 0.5 each)
    m1  = []
    m2 = []
    for i in range(0, k+1):
        m1.append(p1[i][0])
        m2.append(p1[i][1])
    for i in range(k+1, L):
        m1.append(p1[i][1])
        m2.append(p1[i][0])
    m = randint(0,1)
    if m == 0:
        cm = m1
    else:
        cm = m2

    # Choose crossing-over position k’ uniformly among {0, 1, …, L-1}
    k_prime = randint(0, L-1)
    # Child paternal chromosome = c0…ck’ dk’+1…dL-1 or d0…dk’ ck’+1…cL-1 (with probability 0.5 each)
    f1 = []
    f2 = []
    for i in range(0, k_prime+1):
        f1.append(p2[i][0])
        f2.append(p2[i][1])
    for i in range(k_prime+1, L):
        f1.append(p2[i][1])
        f2.append(p2[i][0])
    f = randint(0,1)
    if f == 0:
        cf = f1
    else:
        cf = f2
    child = []
    for i in range(0, L):
        list = []
        list.append(str(cm[i]))
        list.append(str(cf[i]))
        for l in list:
            l = str(l)
        c = "|".join(list)
        child.append(c)
    child_name = "IND"+str(num)

    return pd.Series(child, name=child_name), child_name

def probability_alt_allele_exctinction(population, num_snps):
    """
        probability that a given alternate allele present at frequency 1/2N in generation 0 
        becomes extinct in generation 1. This is simply done by seeing what fraction of the 
        10,000 SNPs went to extinction in generation 1
    """
    num_extinct = 0
    for i in range(0,num_snps-1):
        row = population.iloc[i].values.tolist()
        extinct = is_extinct(row)
        if extinct:
            num_extinct += 1
            
    return num_extinct/num_snps
        

def is_extinct(snp_row):
    arr = []
    extinct = True
    for ind in snp_row:
        s = ind.split('|') 
        s[0] = int(s[0])
        s[1] = int(s[1])
        i = []
        arr.append(s[0])
        arr.append(s[1])
    for i in arr:
        if i == 1:
            extinct = False
            return extinct
    return extinct

def is_fixated(snp_row):
    arr = []
    fixated = True
    for ind in snp_row:
        s = ind.split('|') 
        s[0] = int(s[0])
        s[1] = int(s[1])
        i = []
        arr.append(s[0])
        arr.append(s[1])
    for i in arr:
        if i == 0:
            fixated = False
            return fixated
    return fixated

def probability_alt_allele_fixation(population, num_snps):
    """
        probability that a given alternate allele present at frequency 1/2N in generation 0 
        becomes fixated in generation 1. This is simply done by seeing what fraction of the 
        10,000 SNPs went to fixation in generation 1
    """
    num_fixated = 0
    for i in range(0,num_snps-1):
        row = population.iloc[i].values.tolist()
        fixated = is_fixated(row)
        if fixated:
            num_fixated += 1
    return num_fixated/num_snps


def plot(fitness_func):
    """
        Simulate the evolution of a population over 20 generations.
        Assume that all alleles are selectively neutral.
        Plot the alternate allele frequency of the first 100 SNPs.
    """
    ## Initializations
    num_snps = 100
    num_generations = 20
    gen_0 = pd.read_csv(FILENAME, sep='\t')
    gen_0 = gen_0.drop(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL'], axis=1)

    # Initialize Plot
    #   numpy.zeros returns a new array of given shape and type, filled with zeros.
    arr = np.zeros(shape=(num_generations, num_snps))
    
    # advance evolution and track snps per generation
    parent_gen = gen_0
    for i in range(num_generations):
        arr[i] = count_snps(parent_gen, num_snps)
        parent_gen = evolve(NUM_INDIVIDUALS, num_snps, parent_gen, fitness_func)
    
    print("Plotting...")
    # create a matplotlib plot
    figure, axes = plt.subplots()

    # plot each allele frequency
    for i in range(num_snps):
        axes.plot(range(num_generations), arr[:, i], color=plt.cm.plasma(i/num_snps))
    
    # create scalarMappabe for color range
    norm = colors.Normalize(vmin=0, vmax=num_snps)
    scalar_mappable = cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
    scalar_mappable.set_array([])
    colorbar = plt.colorbar(scalar_mappable, ax=axes, ticks=range(0, 101, 5))
    colorbar.set_label('Allele Index')

    axes.set_ylabel('Alternate Allele Frequency')
    axes.set_xlabel('Generation')
    axes.set_title("Tess Gompper's Evolution of a Population")

    plt.show()

def extinction_fixation_plot(fitness_func):
    num_generations = 1000
    num_snps = 5000
    gen_0 = pd.read_csv(FILENAME, sep='\t')
    gen_0 = gen_0.drop(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL'], axis=1)

    fixations = np.zeros(shape=(num_generations))
    extinctions = np.zeros(shape=(num_generations))
    generations = [0]*num_generations
    for i in range(num_generations):
        generations[i] = i

    # for each generation plot extinciton probability and fixation probability
    parent_gen = gen_0
    for i in range(num_generations):
        if i%100 == 0:
            print(i)
        extinctions[i] = probability_alt_allele_exctinction(parent_gen, num_snps)
        fixations[i] = probability_alt_allele_fixation(parent_gen, num_snps)
        parent_gen = evolve(NUM_INDIVIDUALS, num_snps, parent_gen, fitness_func)
    
    # plot
    plt.figure(figsize=(10, 5))
    plt.plot(generations, extinctions, label='Extinction')
    plt.plot(generations, fixations, label='Fixation')
    plt.title('Probability Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Probability')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.show()
    
def count_snps(population, n):
    snp_arr = [0]*n
    for i in range(0,n-1):
        snp_row = population.iloc[i].values.tolist()
        arr = []
        
        for individual in snp_row:
            s = individual.split('|') 
            s[0] = int(s[0])
            s[1] = int(s[1])
            ind = []
            ind.append(s[0])
            ind.append(s[1])
            arr.append(ind)
        
        for individual in arr:
            if individual[0] == 1 or individual[1] == 1:
                snp_arr[i] += 1/(200)
    return snp_arr        

def snp_42_ext_fixt(fitness_func):
    """
        estimate the probabilities of extinction and fixation of the alternate allele at SNP42, after 100 (not 
        1000) generations. To estimate these probabilities, you will need to run your simulation 
        multiple times (I recommend 1000 times). To make this easier, reduce the number of 
        SNPs to only 100
    """
    num_generations = 100
    num_snps = 100

    pop = population_evolution(NUM_INDIVIDUALS, num_snps, num_generations, fitness_func)
    row = pop.iloc[42].values.tolist()
    extinct = is_extinct(row)
    if extinct:
        ext = 1
    else:
        ext=0
    fixated = is_fixated(row)
    if fixated:
        fixt = 1
    else:
        fixt=0
        #num_extinct += 1
    # ext = probability_alt_allele_exctinction(pop, num_snps)
    # fixt = probability_alt_allele_fixation(pop, num_snps)

    return ext, fixt
    



def snp42_fitness(population, p1=None):
    """
        Return most fit individual.
        Being heterozygous at SNP42 confers a fitness of 1.5 (compared to a fitness of 1 for other alleles)
        and being homozygous for the alternate allele confers a fitness of 2. 
    """
    # isolate row for SNP42
    snp_42 = population.iloc[42].values.tolist()
    population_fitness = 0
    parent_fitness = [0]*NUM_INDIVIDUALS
    i=0
    for individual in snp_42:
        s = individual.split('|') 
        s[0] = int(s[0])
        s[1] = int(s[1])
        ind = []
        ind.append(s[0])
        ind.append(s[1])
        if ind[0]==1 and ind[1] == 1:
            population_fitness+= 2 #0.8 #2
            parent_fitness[i] = 2 #0.8 #2
        elif (ind[0]==1 and ind[1] == 0) or (ind[0]==0 and ind[1] == 1):
            population_fitness+= 1.5 #0.9 #1.5
            parent_fitness[i] = 1.5 #0.9 #1.5
        else:
            population_fitness+=1
            parent_fitness[i] = 1
        i+=1
    if not p1:
        p1_probability_table = [0]*100
        for i in range(100):
            p1_probability_table[i] = parent_fitness[i]/(population_fitness)
        individual = np.random.choice(100, p=p1_probability_table)
    if p1:
        p2_probability_table = [0]*100
        for i in range(100):
            p2_probability_table[i] = parent_fitness[i]/(population_fitness-parent_fitness[p1])
        p2_probability_table[p1]=0
        individual = np.random.choice(100, p=p2_probability_table)
    return individual

def snp42_deleterious_fitness(population, p1=None):
    """
        Return most fit individual.
        Being heterozygous at SNP42 confers a fitness of 1.5 (compared to a fitness of 1 for other alleles)
        and being homozygous for the alternate allele confers a fitness of 2. 
    """
    # isolate row for SNP42
    snp_42 = population.iloc[42].values.tolist()
    population_fitness = 0
    parent_fitness = [0]*NUM_INDIVIDUALS
    i=0
    for individual in snp_42:
        s = individual.split('|') 
        s[0] = int(s[0])
        s[1] = int(s[1])
        ind = []
        ind.append(s[0])
        ind.append(s[1])
        if ind[0]==1 and ind[1] == 1:
            population_fitness+= 0.8
            parent_fitness[i] = 0.8
        elif (ind[0]==1 and ind[1] == 0) or (ind[0]==0 and ind[1] == 1):
            population_fitness+= 0.9 
            parent_fitness[i] = 0.9 
        else:
            population_fitness+=1
            parent_fitness[i] = 1
        i+=1
    if not p1:
        p1_probability_table = [0]*100
        for i in range(100):
            p1_probability_table[i] = parent_fitness[i]/(population_fitness)
        individual = np.random.choice(100, p=p1_probability_table)
    if p1:
        p2_probability_table = [0]*100
        for i in range(100):
            p2_probability_table[i] = parent_fitness[i]/(population_fitness-parent_fitness[p1])
        p2_probability_table[p1]=0
        individual = np.random.choice(100, p=p2_probability_table)
    return individual

def neutral_fitness(population, p1=None):
    """
        Return most fit individual
        fitness is 1 for all
        => return random individual
    """
    N = NUM_INDIVIDUALS
    individual = randrange(N-1)
    if p1:
        while p1 == individual:
            individual=randrange(N-1)
    return individual


if __name__ == "__main__":
    print("running...")
    fitness_func = neutral_fitness
    
    #For question 1c
    # pop = population_evolution(NUM_INDIVIDUALS, NUM_SNPS, 1, neutral_fitness)
    # percent_extinct = probability_alt_allele_exctinction(pop)
    # print(percent_extinct)
    
    # For question 1d
    # plot(fitness_func)

    # For question 1e
    # extinction_fixation_plot(fitness_func)

    # For question 1f
    percent_extinct=0
    num_extinct=0
    for i in range(300):
        if i%50==0:
            print(i)
        pop = population_evolution(NUM_INDIVIDUALS, NUM_SNPS, 1, snp42_fitness)
        row = pop.iloc[42].values.tolist()
        extinct = is_extinct(row)
        if extinct:
            num_extinct += 1
        # percent_extinct += probability_alt_allele_exctinction(pop)..
    percent_extinct=num_extinct/300
    print("f: P(extinction)= " + str(percent_extinct))

    # For question 1g
    ext_sum = 0
    fixt_sum = 0
    for i in range(300): 
        if i%50==0:
            print(i)
        ext, fixt = snp_42_ext_fixt(snp42_fitness)
        ext_sum += ext
        fixt_sum += fixt
    ext_percent = ext_sum/300
    fixt_percent = fixt_sum/300
    print("g: P(fixation)= " + str(fixt_percent))
    print("g: P(extinction)= " + str(ext_percent))

    # # For question 1f
    ext_sum = 0
    fixt_sum = 0
    for i in range(300): 
        if i%50==0:
            print(i)
        ext, fixt = snp_42_ext_fixt(snp42_deleterious_fitness)
        ext_sum += ext
        fixt_sum += fixt
    ext_percent = ext_sum/300
    fixt_percent = fixt_sum/300
    print("f: P(fixation)= " + str(fixt_percent))
    print("f: P(extinction)= " + str(ext_percent))
