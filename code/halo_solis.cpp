// INTRODUCTION

/**
 *  [CASTELLANO]
 * 
 * Practica 4 - Propuesta de Metaheurística: Algoritmo Halo
 * Asignatura: Metaheuristicas
 * Autor: Valentino Lugli (Github: @RhinoBlindado)
 * Fecha: Mayo, Junio 2021
 */

/**
 *  [ENGLISH]
 *
 * Practice 4 - Metaheuristic Proposal: Halo Algorithm
 * Course: Metaheuristics
 * Author: Valentino Lugli (Github: @RhinoBlindado)
 * Fecha: May, June 2021
 */

// LIBRARIES
//  I/O.
#include <iostream>
//  Vectors.
#include <vector>
//  Lists
#include <list>
//  Using the library provided by the professor.
extern "C" 
{
    #include "cec17.h"
}
//  Take the time between executions.
#include <ctime>
//  Do the random shuffle.
#include <algorithm> 
//  Get max size of data types.
#include <limits>
//  Perform math calculations.
#include <cmath>
//  Use the Pair class
#include <utility>
// Randomness
#include <random>

using namespace std;


/*
 *  GLOBAL VARIABLES
 */
const double generationalCrossChance = 0.7,
             mutationChance = 0.1;

const int seed = 19951611;

std::mt19937 predicts(seed);

class cromosome
{
    private:
    vector<double> genes;
    float fitness;
    bool hasChanges;

    public:

    //  Constructors
    cromosome()
    {
        this->fitness = -1;
        this->hasChanges = true;
    }

    cromosome(vector<double> arg_genes, float arg_fitness, bool arg_hasChanges)
    {
        this->genes = arg_genes;
        this->fitness = arg_fitness;
        this->hasChanges = arg_hasChanges;
    }

    // Destructor
    ~cromosome()
    {
        this->genes.clear();
    }

    // Getters
    vector<double> getGenes() const
    {
        return this->genes;
    }

    double gene(int arg_index) const
    {
        return this->genes[arg_index];
    }

    int getGeneSize() const
    {
        return this->genes.size();
    }

    float getFitness() const
    {
        return this->fitness;
    }

    bool getChanges() const
    {
        return this->hasChanges;
    }


    // Setters
    void initGenes(int arg_geneSize)
    {
        this->genes.clear();
        this->genes.resize(arg_geneSize);
    }

    void setGenes(vector<double> arg_genes)
    {
        this->fitness = -1;
        this->hasChanges = true;

        this->genes.clear();
        this->genes = arg_genes;
    }

    void setFitness(float arg_fitness)
    {
        this->fitness = arg_fitness;
    }

    void setChanges(bool arg_hasChanges)
    {
        this->hasChanges = arg_hasChanges;
    }

    void changeGene(int arg_index, double arg_gene)
    {
        this->fitness= -1;
        this->hasChanges = true;

        this->genes[arg_index] = arg_gene;
    }

    void copy(cromosome arg_crom)
    {
        this->fitness = arg_crom.fitness;
        this->hasChanges = arg_crom.hasChanges;
        this->genes = arg_crom.genes;
    }

    // Auxiliar Functions
    void printContents() const
    {
        cout<<"---/Cromosome/---"<<endl;
        cout<<"\tGenes:  [";
        for (int i = 0; i < genes.size(); i++)
        {
            cout<<"i="<<i<<" ["<<genes[i]<<"] ";
        }
        cout<<"]"<<endl;

        cout<<"\tFitness: "<<fitness<<endl;
        cout<<"\tHas changes? "<<hasChanges<<endl;
        cout<<"---/Cromosome End/---"<<endl;
    }

    void printGenes() const
    {
        cout<<"Genes: [ ";
        for (int i = 0; i < genes.size(); i++)
        {
           cout<<"["<<i<<"]="<<genes[i]<<"\t";
        }
        cout<<"]"<<endl;
        
    }
};

struct pop
{
    vector<cromosome> genePool;

    pair<int, double> bestIndividual;
    pair<int, double> worstIndividual;

    double averageFitness;
};

struct bestCromosome
{
    int pos;
    double fitness;
    bool isAlive;
    cromosome data;
};

double inRange(double eval)
{
    return max(-100.0, min(100.0, eval));
}

int evalPop(pop &population)
{
    int popSize = population.genePool.size(),
        counter = 0;
 
    double  actFitness, 
            accumFitness = 0.0,
            bestFitness = numeric_limits<double>::max(),
            worstFitness = numeric_limits<double>::min();

    vector<double> actGenes;

    for (int i = 0; i < popSize; i++)
    {
        if (population.genePool[i].getChanges())
        {
            actGenes = population.genePool[i].getGenes();
            population.genePool[i].setFitness(cec17_fitness(&actGenes[0]));
            population.genePool[i].setChanges(false);
            counter++;
        }

        actFitness = population.genePool[i].getFitness();
        accumFitness += actFitness;

        if (actFitness < bestFitness)
        {
            population.bestIndividual.first = i;
            population.bestIndividual.second = actFitness;
            bestFitness = actFitness;
        }

        if (actFitness > worstFitness)
        {
            population.worstIndividual.first = i;
            population.worstIndividual.second = actFitness;
            worstFitness = actFitness;
        }
    }

    population.averageFitness = accumFitness / popSize;

    return counter;
}

/*
 * SELECTION OPERATOR
 */
void popSelectionGen(const pop &population, vector<cromosome> &parentPop, bestCromosome &bestCrom)
{
    int fighterA,
        fighterB,
        popSize = population.genePool.size(),
        geneSize = population.genePool[0].getGeneSize();    

    std::uniform_int_distribution<> rng(0, popSize-1);

    parentPop = population.genePool;

    bestCrom.isAlive = false;

    for (int i = 0; i < popSize; i++)
    {
        fighterA = rng(predicts);
        fighterB = rng(predicts);
        if (population.genePool[fighterA].getFitness() < population.genePool[fighterB].getFitness())
        {
            parentPop[i] = population.genePool[fighterA];

            if (fighterA == bestCrom.pos)
            {
                bestCrom.pos = i;
                bestCrom.isAlive = true;
            }
        }
        else
        {
            parentPop[i] = population.genePool[fighterB];
            if (fighterB == bestCrom.pos)
            {
                bestCrom.pos = i;
                bestCrom.isAlive = true;
            }
        }
    }
}

/*
 *  CROSSOVER OPERATORS
 */

void BLX_Alpha(vector<cromosome> &population, bestCromosome &bestCrom)
{
    int crossNum = (int)ceil(population.size() * generationalCrossChance), geneSize = population[0].getGeneSize();

    vector<double> child1(geneSize), child2(geneSize);

    double cMin, cMax, I, lowBound, highBound, alpha = 0.5;

    if (bestCrom.pos < crossNum)
    {
        bestCrom.isAlive = false;
    }

 //   cout << "Starting main loop" << endl;
 //   cout << "CrossNum=" << crossNum << " Genesize=" << geneSize << endl;
    for (int i = 0; i < crossNum; i += 2)
    {
        for (int j = 0; j < geneSize; j++)
        {
          //  cout << "Calculating BLX" << endl;
            cMin = min(population[i].gene(j), population[i + 1].gene(j));
            cMax = max(population[i].gene(j), population[i + 1].gene(j));
            I = cMax - cMin;
            lowBound = inRange(cMin - I * alpha);
            highBound = inRange(cMax + I * alpha);
            std::uniform_real_distribution<> blx(lowBound, highBound);
            child1[j] = blx(predicts);
            child2[j] = blx(predicts);
          //  cout << "BLX DONE" << endl;
        }
        population[i].setGenes(child2);
        population[i + 1].setGenes(child1);
    }
}


void arithmetic(vector<cromosome> &population, bestCromosome &bestCrom)
{
    int crossNum = (int)ceil(population.size() * generationalCrossChance), geneSize = population[0].getGeneSize();

    vector<double> child1(geneSize), child2(geneSize);
    
    double lambda;
    
    std::uniform_real_distribution<> lambdaGen(0.25, 0.75);

    if (bestCrom.pos < crossNum)
    {
        bestCrom.isAlive = false;
    }

    for (int i = 0; i < crossNum; i += 2)
    {
        for (int j = 0; j < geneSize; j++)
        {
            lambda = lambdaGen(predicts);

            child1[j] = lambda * population[i].gene(j) + (1 - lambda) * population[i + 1].gene(j);
            child2[j] = lambda * population[i + 1].gene(j) + (1 - lambda) * population[i].gene(j);
        }
        population[i].setGenes(child2);
        population[i + 1].setGenes(child1);
    }
}


bool compareFitness(const pair<double, vector<double>> &a, const pair<double, vector<double>> &b)
{
    return a.first < b.first;
}

int lineal(vector<cromosome> &population, bestCromosome &bestCrom)
{
    int crossNum = (int)ceil(population.size() * generationalCrossChance), geneSize = population[0].getGeneSize(), localEvals = 0;

    double actCalc;

    vector<double> child1(geneSize), child2(geneSize), child3(geneSize);
    list<pair<double, vector<double>>> floodPod;

    if (bestCrom.pos < crossNum)
    {
        bestCrom.isAlive = false;
    }

    for (int i = 0; i < crossNum; i+=2)
    {
        floodPod.clear();
        for (int j = 0; j < geneSize; j++)
        {
            actCalc = 0.5 * population[i].gene(j) + 0.5 * population[i + 1].gene(j);
            child1[j] = inRange(actCalc);

            actCalc = 1.5 * population[i].gene(j) - 0.5 * population[i + 1].gene(j);
            child2[j] = inRange(actCalc);

            actCalc = -0.5 * population[i].gene(j) + 1.5 * population[i + 1].gene(j); 
            child3[j] = inRange(actCalc);
        }

        floodPod.push_back(pair<double, vector<double>>(cec17_fitness(&child1[0]), child1));
        floodPod.push_back(pair<double, vector<double>>(cec17_fitness(&child2[0]), child2));
        floodPod.push_back(pair<double, vector<double>>(cec17_fitness(&child3[0]), child3));

        localEvals += 3;

        floodPod.sort(compareFitness);

        population[i].setGenes(floodPod.front().second);
        population[i].setFitness(floodPod.front().first);
        population[i].setChanges(false);
        floodPod.pop_front();

        population[i + 1].setGenes(floodPod.front().second);
        population[i + 1].setFitness(floodPod.front().first);
        population[i + 1].setChanges(false);
        floodPod.pop_front();
    }

    return localEvals;
}

/*
 * MUTATION OPERATORS
 */
double delta_nonUniform(double t, double y, double r, double gMax)
{
    return y * (1 - pow(r, pow(1 - (t / gMax), 2)));
}

void nonUniformMutation(vector<cromosome> &population, bestCromosome &bestCrom, int actGen, int maxGen)
{
    int mutationSize = mutationChance * population.size(), ii, jj;
    double actGene;
    std::uniform_int_distribution<int> tau(0, 1), startCromosome(0, population.size() - 1 - mutationSize), mutateGene(0, population[0].getGeneSize() - 1);
    std::uniform_real_distribution<double> r(0, 1);

    ii = startCromosome(predicts);

    if (bestCrom.pos >= ii && bestCrom.pos <= ii+mutationSize)
    {
        bestCrom.isAlive = false;
    }

    
    for (int i = 0; i < mutationSize; i++)
    {
        jj = mutateGene(predicts);

        actGene = population[i + ii].gene(jj);

        if (tau(predicts))
        {
            actGene = actGene - delta_nonUniform(actGen, actGene - 100.0, r(predicts), maxGen);
        }
        else
        {
            actGene = actGene + delta_nonUniform(actGen, 100.0 - actGene, r(predicts), maxGen);
        }

        population[i + ii].changeGene(jj, inRange(actGene));

    }

}

void muhlenbeinMutation(vector<cromosome> &population, bestCromosome &bestCrom)
{
    int mutationSize = mutationChance * population.size(), ii, jj;
    double actGene, gamma;
    std::uniform_int_distribution<int> plusOrMinus(0, 1), startCromosome(0, population.size() - 1 - mutationSize), mutateGene(0, population[0].getGeneSize() - 1), alpha(0, 15);
    std::uniform_real_distribution<double> rang(-10.0, 10.0);

    ii = startCromosome(predicts);

    if (bestCrom.pos >= ii && bestCrom.pos <= ii+mutationSize)
    {
        bestCrom.isAlive = false;
    }

    for (int i = 0; i < mutationSize; i++)
    {
        jj = mutateGene(predicts);

        actGene = population[i + ii].gene(jj);

        gamma = 0.0;

        for (int k = 0; k < 15; k++)
        {
            if (alpha(predicts) == 1)
            {
                gamma += k * pow(2, -k);
            }
        }

        if (plusOrMinus(predicts))
        {
            actGene += rang(predicts) * gamma;
        }
        else
        {
            actGene -= rang(predicts) * gamma;
        }
        population[i + ii].changeGene(jj, inRange(actGene));
    }
} 

void randomMutation(vector<cromosome> &population, bestCromosome &bestCrom)
{
    int mutationSize = mutationChance * population.size(), ii, jj;
    double actGene;
    std::uniform_int_distribution<int> startCromosome(0, population.size() - 1 - mutationSize), mutateGene(0, population[0].getGeneSize() - 1);
    std::uniform_real_distribution<double> floodSpore(-100.0, 100.0);

    ii = startCromosome(predicts);

    if (bestCrom.pos >= ii && bestCrom.pos <= ii+mutationSize)
    {
        bestCrom.isAlive = false;
    }

    for (int i = 0; i < mutationSize; i++)
    {
        jj = mutateGene(predicts);

        actGene = floodSpore(predicts);

        population[i + ii].changeGene(jj, actGene);
    }
}

/*
 * REPLACEMENT OPERATOR
 */
void replacePopGen(pop &childPop, vector<cromosome> parentPop, bestCromosome bestCrom)
{
    for (int i = 0; i < parentPop.size(); i++)
    {
        childPop.genePool[i] = parentPop[i];        
    }

    if (!bestCrom.isAlive)
    {
        childPop.genePool[0] = bestCrom.data;
    }

}

int AGG(pop &population, int dim, int popType)
{
    bool localOptimaFound = false;

    int evaluations = 0,
        localOptimaCounter = 0,
        generations = 0;

    evaluations += evalPop(population);

    double prevFitness = population.bestIndividual.second;

    vector<cromosome> parentPop;
    bestCromosome bestCrom;

    parentPop.resize(population.genePool.size());

    while (evaluations < 250 * dim)
    {
      //  cout << "POP INFO: " << population.bestIndividual.first << endl;
        bestCrom.pos = population.bestIndividual.first;
        bestCrom.fitness = population.bestIndividual.second;
        bestCrom.data = population.genePool[bestCrom.pos];
        bestCrom.isAlive = true;

   /*     cout << "Iteracion " << evaluations << endl;
        for (int i = 0; i < population.genePool.size(); i++)
        {
            cout << "[" << i << "]";
            population.genePool[i].printGenes();
        }
        cout << "El mejor es " << bestCrom.pos <<"("<<bestCrom.fitness<<")"<<endl;


            cout << "SELECTION OPERATOR" << endl;*/
        popSelectionGen(population, parentPop, bestCrom);

        switch (popType)
        {
            // Humanity
            case 0:
                BLX_Alpha(parentPop, bestCrom);
                nonUniformMutation(parentPop, bestCrom, evaluations, 250*dim);
            break;

            // The Covenant
            case 1:
                arithmetic(parentPop, bestCrom);
                muhlenbeinMutation(parentPop, bestCrom);
            break;

            // The Flood
            case 2:
                evaluations += lineal(parentPop, bestCrom);
                randomMutation(parentPop, bestCrom);
            break;
        }
  //      if(dim == 30)
         //   cout << "REPLACEMENT OPERATOR" << endl;
        replacePopGen(population, parentPop, bestCrom);

        prevFitness = population.bestIndividual.second;

        evaluations += evalPop(population);

     //   cout << "MEJOR GEN ES " << population.bestIndividual.first << " CON FITNESS " << population.bestIndividual.second << endl;

        if (prevFitness <= population.bestIndividual.second)
        {
            localOptimaCounter++;
      //      cout << "Estancado... " << localOptimaCounter << " veces" << endl;
            if (localOptimaCounter >= (250 * dim) / 2)
            {
                break;
            }
        }
        else
        {
     //       cout << "Mejora" << endl;
            localOptimaCounter = 0;
        }
    }

   return evaluations;
}


void initPop(pop &population, int popSize, int geneDim)
{
    std::uniform_real_distribution<double> theOracle(-100.0, 100.0);

    population.genePool.resize(popSize);
    
    for (int i = 0; i < popSize; i++)
    {
        population.genePool[i].initGenes(geneDim);

        for (int j = 0; j < geneDim; j++)
        {
            population.genePool[i].changeGene(j, theOracle(predicts));
        }
    }
}

vector<double> regenerate(const vector<double> &parentA, const vector<double> &parentB)
{
    double a, b;
    vector<double> newChild(parentA.size());

    for (int i = 0; i < newChild.size(); i++)
    {
        a = parentA[i];
        b = parentB[i];
        std::uniform_real_distribution<double> genesis(min(a, b), max(a, b));
        newChild[i] = genesis(predicts);
    }

    return newChild;
}

void popAdvantage(pop &mainPop, pop &primaryPop, pop &secondaryPop)
{
    std::uniform_int_distribution<int> lengthPrimary((int)mainPop.genePool.size() * 0.25, (int)mainPop.genePool.size() * 0.5),
                                       lengthSecondary(0, (int)mainPop.genePool.size() * 0.5),
                                       start(0, mainPop.genePool.size() - 1);

    int sP, fP, sS, fS;

    sP = start(predicts);
    fP = lengthPrimary(predicts);

    sS = start(predicts);
    fS = lengthSecondary(predicts);

    for (int i = sP; i < fP; i++)
    {
        if (primaryPop.genePool[i].getFitness() > mainPop.genePool[i].getFitness())
        {
            primaryPop.genePool[i].setGenes(regenerate(primaryPop.genePool[primaryPop.bestIndividual.first].getGenes(), primaryPop.genePool[start(predicts)].getGenes()));
        }
        else
        {
            mainPop.genePool[i].copy(primaryPop.genePool[i]);
        }
    }

    for (int i = sS; i < fS; i++)
    {
        if (secondaryPop.genePool[i].getFitness() > mainPop.genePool[i].getFitness())
        {
            secondaryPop.genePool[i].setGenes(regenerate(secondaryPop.genePool[secondaryPop.bestIndividual.first].getGenes(), secondaryPop.genePool[start(predicts)].getGenes()));
        }
        else
        {
            mainPop.genePool[i].copy(secondaryPop.genePool[i]);
        }
    }
    
}

// SOLIS WETS FUNCTIONS

void clip(vector<double> &sol, int lower, int upper) {
  for (auto &val : sol) {
    if (val < lower) {
      val = lower;
    }
    else if (val > upper) {
      val = upper;
    }
  }
}

void increm_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = 0.2*bias[i]+0.4*(dif[i]+bias[i]);
  }
}

void decrement_bias(vector<double> &bias, vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = bias[i]-0.4*(dif[i]+bias[i]);
  }
}

/**
 * Aplica el Solis Wets
 *
 * @param  sol solucion a mejorar.
 * @param fitness fitness de la solución.
 */
template <class Random>
void soliswets(vector<double> &sol, double &fitness, double delta, int maxevals, int lower, int upper, Random &random) {
    const size_t dim = sol.size();
    vector<double> bias (dim), dif (dim), newsol (dim);
    double newfit;
    size_t i;

    int evals = 0;
    int num_success = 0;
    int num_failed = 0;

    while (evals < maxevals) {
        std::uniform_real_distribution<double> distribution(0.0, delta);

        for (i = 0; i < dim; i++) {
        dif[i] = distribution(random);
        newsol[i] = sol[i] + dif[i] + bias[i];
        }

        clip(newsol, lower, upper);
        newfit = cec17_fitness(&newsol[0]);
        evals += 1;

        if (newfit < fitness) {
        sol = newsol;
        fitness = newfit;
        increm_bias(bias, dif);
        num_success += 1;
        num_failed = 0;
        }
        else if (evals < maxevals) {

        for (i = 0; i < dim; i++) {
            newsol[i] = sol[i] - dif[i] - bias[i];
        }

        clip(newsol, lower, upper);
        newfit = cec17_fitness(&newsol[0]);
        evals += 1;

        if (newfit < fitness) {
            sol = newsol;
            fitness = newfit;
            decrement_bias(bias, dif);
            num_success += 1;
            num_failed = 0;
        }
        else {
            for (i = 0; i < dim; i++) {
            bias[i] /= 2;
            }

            num_success = 0;
            num_failed += 1;
        }
        }

        if (num_success >= 5) {
        num_success = 0;
        delta *= 2;
        }
        else if (num_failed >= 3) {
        num_failed = 0;
        delta /= 2;
        }
    }
}


void halo(int dim)
{
    int evaluations = 0,
        years = 0,
        popSize = 20;

    double bestOverallFitness, auxFitness;

    vector<double> auxGenes;

    pop humanity, 
        theCovenant, 
        theFlood;

    initPop(humanity, popSize, dim);
    initPop(theCovenant, popSize, dim);
    initPop(theFlood, popSize, dim);
    
    evalPop(humanity);
    evalPop(theCovenant);
    evalPop(theFlood);
    
    do
    {
        // War Phase
        if (years == 2)
        {
            //cout << "WAR" << endl;
            years = 0;
            if (bestOverallFitness == humanity.bestIndividual.second)
            {
                if (theCovenant.averageFitness < theFlood.averageFitness)
                {
                    popAdvantage(humanity, theFlood, theCovenant);
                }
                else
                {
                    popAdvantage(humanity, theCovenant, theFlood);
                }
            }
            else if (bestOverallFitness == theCovenant.bestIndividual.second)
            {
                if (humanity.averageFitness < theFlood.averageFitness)
                {
                    popAdvantage(theCovenant, theFlood, humanity);
                }
                else
                {
                    popAdvantage(theCovenant, humanity, theFlood);
                }
            }
            else
            {
                if (humanity.averageFitness < theCovenant.averageFitness)
                {
                    popAdvantage(theFlood, theCovenant, humanity);
                }
                else
                {
                    popAdvantage(theFlood, humanity, theCovenant);
                }
            }
           // cout << "WAR END" << endl;
            auxGenes = humanity.genePool[humanity.bestIndividual.first].getGenes();
            auxFitness = humanity.bestIndividual.second;
            
            soliswets(auxGenes, auxFitness, 0.2, 150 * dim , -100, 100, predicts);

            if (auxFitness < humanity.bestIndividual.second)
            {
                humanity.genePool[humanity.bestIndividual.first].setGenes(auxGenes);
                humanity.genePool[humanity.bestIndividual.first].setFitness(auxFitness);
                humanity.genePool[humanity.bestIndividual.first].setChanges(false);
                humanity.bestIndividual.second = auxFitness;
            }

            auxGenes = theCovenant.genePool[theCovenant.bestIndividual.first].getGenes();
            auxFitness = theCovenant.bestIndividual.second;
            
            soliswets(auxGenes, auxFitness, 0.2, 150 * dim, -100, 100, predicts);

            if (auxFitness < theCovenant.bestIndividual.second)
            {
                theCovenant.genePool[theCovenant.bestIndividual.first].setGenes(auxGenes);
                theCovenant.genePool[theCovenant.bestIndividual.first].setFitness(auxFitness);
                theCovenant.genePool[theCovenant.bestIndividual.first].setChanges(false);
                theCovenant.bestIndividual.second = auxFitness;
            }

            auxGenes = theFlood.genePool[theFlood.bestIndividual.first].getGenes();
            auxFitness = theFlood.bestIndividual.second;
            
            soliswets(auxGenes, auxFitness, 0.2, 150 * dim, -100, 100, predicts);

            if (auxFitness < theFlood.bestIndividual.second)
            {
                theFlood.genePool[theFlood.bestIndividual.first].setGenes(auxGenes);
                theFlood.genePool[theFlood.bestIndividual.first].setFitness(auxFitness);
                theFlood.genePool[theFlood.bestIndividual.first].setChanges(false);
                theFlood.bestIndividual.second = auxFitness;
            }
        }
        years++;

        // Peace Phase
        evaluations += AGG(humanity, dim, 0);
        evaluations += AGG(theCovenant, dim, 1);
        evaluations += AGG(theFlood, dim, 2);

        bestOverallFitness = min({humanity.bestIndividual.second, theCovenant.bestIndividual.second, theFlood.bestIndividual.second});
    } while (evaluations < 10000 * dim);
    cout << "\t Fitness = " << bestOverallFitness << endl;
}

/**
 * @brief Main function
 * @return Execution status
 */
int main()
{
    // Initialice parameters
    //  - Dimentions to be evaluated.
    vector<int> dimension = {10, 30};

    // Start the main loop.
    cout << "-- Halo x Solis Wets Memetic Algorithm --" << endl;

    for (int i = 0; i < dimension.size(); i++)
    {
        cout << ">>Dimension: " << dimension[i] << endl;
        for (int funcid = 1; funcid <= 30; funcid++)
        {
            cout << ">>>Function " << funcid << endl;
            for (int k = 0; k < 10; k++)
            {
                cout << "\t Execution [" << k + 1 << "/10]...\n";
                cec17_init("HaloSolis", funcid, dimension[i]);
                halo(dimension[i]);
            }
        }
    }

    return 0;
}
