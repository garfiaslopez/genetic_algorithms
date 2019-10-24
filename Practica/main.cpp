// compile & run:     g++ main.cpp -o main && ./main > out.txt
// © Jose Garfias Lopez
// Genetic  Algorithm

#include <iostream>
#include <cstdio>
#include <vector>
#include <math.h>

using namespace std; 

#define POPULATION_SIZE 100
#define SIZE_OF_CHROMOSOME 390
#define SIZE_OF_CHUNK 13
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
    void calc_fitness(); 
    void uniform_mutate(); 
};
  
Individual::Individual(string chromosome) {
    this->chromosome = chromosome; 
    this->calc_fitness();  
}; 

vector<Individual> Individual::mate(Individual partner) {
    vector<Individual> childs;
    string parentA_A = "";
    string parentA_B = "";
    string parentB_A = "";
    string parentB_B = "";
    int point_of_cross = random_num(1, SIZE_OF_CHROMOSOME - 1);

    for (int i=0; i <= point_of_cross; i++) {
        parentA_A += this->chromosome[i];
        parentB_A += partner.chromosome[i];
    }
    for (int i=point_of_cross+1; i < SIZE_OF_CHROMOSOME; i++) {
        parentA_B += this->chromosome[i];
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
    float aj = -30;
    float bj = 30;
    return aj + (float)binary_to_decimal(chromosome) * ((bj - aj) / ((pow(2,13) - 1)));
}

float fitness_func(string chromosome) {
    // CALCULE X PER SUBCHROMOSOME IN CHROMOSOME
    vector<float> x_values;
    int chunk = 0;
    string sub_chromosome = "";
    for (int i=0; i<SIZE_OF_CHROMOSOME; i++) {
        sub_chromosome += chromosome[i];
        if (sub_chromosome.length() == 13) {
            x_values.push_back(chromosome_to_x(sub_chromosome));
            sub_chromosome = "";
        }
    }

    // EVAL FITNESS FUNCTION
    float sum_result = 0;
    for (int i=0; i<x_values.size(); i++) {
        float Xi = x_values[i];
        sum_result += abs(100 * powf((Xi+1)-(powf(Xi, 2)), 2) + powf((Xi-1), 2));
    }
    return sum_result;
}

void Individual::calc_fitness() {
    this->fitness = fitness_func(this->chromosome);
};

void Individual::uniform_mutate() {
    int random_position = random_num(0, SIZE_OF_CHROMOSOME - 1);
    char random_gene = mutated_genes();
    this->chromosome[random_position] = random_gene;
    this->calc_fitness();
};
  
bool operator<(const Individual &ind1, const Individual &ind2) {
    return ind1.fitness > ind2.fitness; 
}

Individual roulete_selection (vector<Individual> population) {
    int n = population.size();
    vector<float> p_slect_array;
    vector<float> val_esp_array;
    float T = 0;
    float sumOfFitness = 0;
    for (int i=0; i<n; i++) {
        sumOfFitness += population[i].fitness;
    }
    for (int i=0; i<n; i++) {
        p_slect_array.push_back(1 - ((population[i].fitness) / sumOfFitness));
        val_esp_array.push_back(p_slect_array[i] * n);
        T += val_esp_array[i];
    }

    float r = random_float_num(0.0, T);
    int individual_i = -1;

    float sum=0;
    for (int i=0; i<n; i++) {
        sum += val_esp_array[i];
        if (sum >= r) {
            individual_i = i;
            break;
        }
    }

    return Individual(population[individual_i].chromosome); 
}

void print_generation(int genration_number, vector<Individual> population, bool best_member) {
    if (best_member) {
        cout<< "Generation: " << genration_number << "\t"; 
        cout<< "String: "<< population[0].chromosome <<"\t"; 
        cout<< "Value: "<< chromosome_to_x(population[0].chromosome) <<"\t"; 
        cout<< "Fitness: "<< population[0].fitness << "\n";
    } else {
        for (int i=0; i<population.size(); i++) {
            cout<< "Generation: " << genration_number << "\t"; 
            cout<< "String: "<< population[i].chromosome <<"\t"; 
            cout<< "Value: "<< chromosome_to_x(population[i].chromosome) <<"\t"; 
            cout<< "Fitness: "<< population[i].fitness << "\n";
        }
    }
    
}
  
int main() {
    srand((unsigned)(time(0)));
    int generation = 0; 
    vector<Individual> population; 
    vector<Individual> new_generation;

    // create initial population
    for(int i=0; i<POPULATION_SIZE; i++) {
        string gnome = create_gnome();
        population.push_back(Individual(gnome));
    }
  
    while(generation <= 3000) {
        new_generation.clear(); // clear the array for new generation

        // SORT THE GENERATION BY FITNESS
        sort(population.begin(), population.end()); 
        print_generation(generation, population, true);

        // ELITISM
        int s = (20*POPULATION_SIZE)/100; 
        for(int i = 0; i<s; i++){
            new_generation.push_back(population[i]);
        }
    
        // SELECTION AND CROSSOVER
        s = ((80*POPULATION_SIZE)/100 ) / 2;
        for(int i = 0; i<s; i++) {
            Individual parent1 = roulete_selection(population);
            Individual parent2 = roulete_selection(population);
            vector<Individual> childs = parent1.mate(parent2);
            new_generation.push_back(childs[0]);
            new_generation.push_back(childs[1]);
        }

        // MUTATION: 1/L
        new_generation[0].uniform_mutate(); // the best

        population = new_generation;
        generation++;
     }
}