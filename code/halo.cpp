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
//  Using the library provided by the professor.
extern "C" {
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

using namespace std;

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

    void changeGene(int arg_index, int arg_gene)
    {
        this->fitness= -1;
        this->hasChanges = true;

        this->genes[arg_index] = arg_gene;
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


void UNSC()
{

}

void theConvenant()
{

}

void theFlood()
{

}

void halo(int dim)
{
    
}


/**
 * @brief Main function
 * @return Execution status
 */
int main()
{
    vector<int> dimension = {10, 30};
    cout << "-- Halo Genetic Algorithm --" << endl;

    for (int i = 0; i < dimension.size(); i++)
    {
        cout << ">>Dimension: " << dimension[i] << endl;
        for (int funcid = 1; funcid <= 30; funcid++)
        {
            cout << ">>>Function " << funcid << endl;
            for (int k = 0; k < 10; k++)
            {
                cout << "\t Execution [" << k + 1 << "/10]...";

                cec17_init("Halo", funcid, dimension[i]);
                halo(dimension[i]);
                cout << "Done" << endl;
            }
        }
    }

    return 0;
}
