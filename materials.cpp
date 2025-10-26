#include "materials.h"

#include <iostream>

using namespace std;

FunctionElement::FunctionElement(double x, double fx) {
    this->x = x;
    this->fx = fx;
}

FunctionElement::FunctionElement() : FunctionElement(0, 0) {
}

double FunctionElement::getX() {
    return this->x;
}

double FunctionElement::getFx() {
    return this->fx;
}

void FunctionElement::multiplyXBy(int k) {
    x *= k;
}

void FunctionElement::multiplyFxBy(int k) {
    fx *= k;
}

FunctionElement::~FunctionElement() {
}

FunctionSpectre::FunctionSpectre() : FunctionSpectre(0) {
}

FunctionSpectre::FunctionSpectre(int spectreLen) {
    this->spectreLen = spectreLen;
    spectre = new FunctionElement[spectreLen];
}

double FunctionSpectre::selectValue(double x) {
    double fx = 0;
    if(spectreLen > 0) {
        int i = 0;
        while(i < spectreLen && spectre[i].getX() < x) i++;
        if(i >= spectreLen) fx = spectre[spectreLen - 1].getFx();
        else if(i == 0 || spectre[i].getX() == x || spectre[i].getX() - spectre[i - 1].getX() == 0) fx = spectre[i].getFx();
        else fx = spectre[i - 1].getFx() + ((spectre[i].getFx() - spectre[i - 1].getFx()) / (spectre[i].getX() - spectre[i - 1].getX())) * (x - spectre[i - 1].getX());
    } else {
        cout << "ERROR: empty spectre" << endl;
    }
    return fx;
}

FunctionSpectre** FunctionSpectre::readMultipleFromFile(int len, std::string filePath, bool invertAndNormalize, double k) {
    FILE* ptr = fopen(filePath.c_str(), "r");
    int spectreLen = 0;
    FunctionSpectre** spects = nullptr;
	if(ptr == NULL) {
        cout << "ERROR: Cannot open: " << filePath << endl;
	}
	else {
        double temp;
        while(fscanf(ptr, "%lf", &temp) != EOF) {
            for(int i = 0; i < len; i++) fscanf(ptr, "%lf", &temp);
            spectreLen++;
        }
        spects = new FunctionSpectre*[len];
        for(int i = 0; i < len; i++) {
            spects[i] = new FunctionSpectre(spectreLen);
        }

        fseek(ptr, 0, SEEK_SET);
        for(int i = 0; i < spectreLen; i++) {
            double x, fx;
            fscanf(ptr, "%lf", &x);
            for(int j = 0; j < len; j++) {
                fscanf(ptr, "%lf", &fx);
                spects[j]->spectre[i] = FunctionElement(x, fx);
            }
        }

        if(invertAndNormalize) {
            for(int i = 0; i < len; i++) {
                double s = 0;
                for(int j = 0; j < spectreLen; j++) {
                    s += spects[i]->spectre[j].getFx();
                }
                if(s != 0) {
                    for(int j = 0; j < spectreLen; j++) {
                        if(j == 0) {
                            spects[i]->spectre[j] = FunctionElement(spects[i]->spectre[j].getFx() / s, spects[i]->spectre[j].getX());
                        } else {
                            spects[i]->spectre[j] = FunctionElement(spects[i]->spectre[j].getFx() / s + spects[i]->spectre[j - 1].getX(), spects[i]->spectre[j].getX());
                        }

                    }
                }
            }
        } else {
            for(int i = 0; i < len; i++) {
                for(int j = 0; j < spectreLen; j++) {
                    spects[i]->spectre[j] = FunctionElement(1000 * spects[i]->spectre[j].getX(), spects[i]->spectre[j].getFx() * k);
                }
            }
        }
	}
	fclose(ptr);
	return spects;
}

FunctionSpectre* FunctionSpectre::readFromFile(std::string filePath, bool invertAndNormalize, double k) {
    auto spects = FunctionSpectre::readMultipleFromFile(1, filePath, invertAndNormalize, k);
    if(spects != nullptr) {
        auto spectre = spects[0];
        spects[0] = nullptr;
        delete[] spects;
        return spectre;
    }

    return nullptr;
}

FunctionSpectre::~FunctionSpectre() {
    delete[] spectre;
}
