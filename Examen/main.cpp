// compile & run:     g++ main.cpp -o main && ./main > out.txt
// © Jose Garfias Lopez
// Genetic  Algorithm

#include <iostream>
#include <cstdio>
#include <vector>
#include <math.h>

using namespace std; 

#define POPULATION_SIZE 50
#define SIZE_OF_CHROMOSOME 13
const string GENES = "01"; 

int random_num(int start, int end) {
    int range = (end - start) + 1; 
    int random_int = start + (rand() % range); 
    return random_int; 
}

float random_float_num(float start, float end) {
    int range = (end - start) + 1; 
    float random_float = start + (rand() % range); 
    return random_float; 
}
  
char mutated_genes() { 
    int r = random_num(0, GENES.size() - 1); 
    return GENES[r]; 
} 
  
string create_gnome() {
    string gnome = "";
    for(int i=0; i < SIZE_OF_CHROMOSOME; i++) {
        gnome += mutated_genes();
    }
    return gnome; 
} 
  
// Class representing individual in population 
class Individual {
public: 
    string chromosome; 
    float fitness; 
    Individual(string chromosome); 
    vector<Individual> mate(Individual parent2); 
    float cal_fitness(); 
};
  
Individual::Individual(string chromosome) {
    this->chromosome = chromosome; 
    fitness = cal_fitness();  
}; 

vector<Individual> Individual::mate(Individual partner) {
    vector<Individual> childs;
    string parentA_A = "";
    string parentA_B = "";
    string parentB_A = "";
    string parentB_B = "";
    int point_of_cross = random_num(1, SIZE_OF_CHROMOSOME - 1);

    for (int i=0; i <= point_of_cross; i++) {
        parentA_A += chromosome[i];
        parentB_A += partner.chromosome[i];
    }
    for (int i=point_of_cross+1; i < SIZE_OF_CHROMOSOME; i++) {
        parentA_B += chromosome[i];
        parentB_B += partner.chromosome[i];
    }

    childs.push_back(Individual(parentA_A + parentB_B));
    childs.push_back(Individual(parentA_B + parentB_A));

    return childs; 
};

int binary_to_decimal(string n) {
    string num = n; 
    int dec_value = 0; 
    int base = 1;  
    int len = num.length(); 
    for (int i = len - 1; i >= 0; i--) {
        if (num[i] == '1') {
            dec_value += base; 
        }
        base = base * 2; 
    }
    return dec_value; 
} 

float chromosome_to_x(string chromosome) {
    float aj = -5;
    float bj = 5;
    return aj + (float)binary_to_decimal(chromosome) * ((bj - aj) / ((pow(2,13) - 1)));
}

float fitness_func(float x) {
    return 0.65 - (0.75/(1 + powf(x, 2))) - (0.65 * x * atanf(1/x));
}

float Individual::cal_fitness() {
    return fitness_func(chromosome_to_x(chromosome));
};
  
bool operator<(const Individual &ind1, const Individual &ind2) {
    return ind1.fitness < ind2.fitness; 
}

Individual roulete_selection (vector<Individual> population) {
    int n = population.size();
    vector<float> p_slect_array;
    vector<float> val_esp_array;
    float T = 0;
    float sumOfFitness = 0;
    for (int i = 0; i<n; i++) {
        sumOfFitness += population[0].fitness;
    }
    for (int i = 0; i<n; i++) {
        p_slect_array[i] = 1 - ((population[0].fitness) / sumOfFitness);
        val_esp_array[i] = p_slect_array[i] * n;
        T += val_esp_array[i];
    }

    float r = random_float_num(0.0, T);
    int individual_i = -1;

    for (int i = 0; i<n; i++) {
        float sum = 0;
        sum += val_esp_array[0];
        if (sum >= r) {
            individual_i = i;
            break;
        }
    }
    
    return Individual(population[individual_i].chromosome); 
}
  
int main() {
    srand((unsigned)(time(0)));
    int generation = 0; 
    vector<Individual> population; 

    // create initial population
    for(int i=0; i<POPULATION_SIZE; i++) {
        string gnome = create_gnome();
        population.push_back(Individual(gnome));
    } 
  
    while(generation < 50) {
        sort(population.begin(), population.end()); 
  
        vector<Individual> new_generation; 
  
        // ELITISM
        int s = (20*POPULATION_SIZE)/100; 
        for(int i = 0; i<s; i++){
            new_generation.push_back(population[i]);
        }
    
        // CROSSOVER
        s = ((80*POPULATION_SIZE)/100 ) / 2;
        for(int i = 0; i<s; i++) {
            int r = random_num(0, POPULATION_SIZE - 1);
            Individual parent1 = population[r];
            r = random_num(0, POPULATION_SIZE - 1);
            Individual parent2 = population[r];
            vector<Individual> childs = parent1.mate(parent2);
            new_generation.push_back(childs[0]);
            new_generation.push_back(childs[1]);
        }

        population = new_generation;

        
        cout<< "Generation: " << generation << "\t"; 
        cout<< "String: "<< population[0].chromosome <<"\t"; 
        cout<< "Value: "<< chromosome_to_x(population[0].chromosome) <<"\t"; 
        cout<< "Fitness: "<< population[0].fitness << "\n"; 
  
        generation++; 
     } 
    
    cout<< "Generation: " << generation << "\t"; 
    cout<< "String: "<< population[0].chromosome <<"\t"; 
    cout<< "Value: "<< chromosome_to_x(population[0].chromosome) <<"\t"; 
    cout<< "Fitness: "<< population[0].fitness << "\n"; 
}