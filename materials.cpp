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

FunctionElement::~FunctionElement() {
}

FunctionSpectre::FunctionSpectre() {
    spectreLen = 0;
    spectre = new FunctionElement[spectreLen];
}

FunctionSpectre::FunctionSpectre(std::string filePath, FunctionSpectreOrigin origin, double k) {
    FILE* ptr = fopen(filePath.c_str(), "r");
    spectreLen = 0;
	if(ptr == NULL) {
        spectre = new FunctionElement[spectreLen];
        cout << "ERROR: Cannot open: " << filePath << endl;
	}
	else {
        double temp, s;
        switch(origin) {
            case FunctionSpectreOrigin::ENERGY: {
                s = 0;
                while (fscanf(ptr, "%lf %lf", &temp, &temp) != EOF) {
                    spectreLen++;
                    s += temp;
                }
                break;
            }
            case FunctionSpectreOrigin::MATERIAL: {
                while (fscanf(ptr, "%lf %lf %lf", &temp, &temp, &temp) != EOF) {
                    spectreLen++;
                }
                break;
            }
        }

        fseek(ptr, 0, SEEK_SET);
        spectre = new FunctionElement[spectreLen];
        int i = 0;
        double tmp1, tmp2;
        switch(origin) {
            case FunctionSpectreOrigin::ENERGY: {
                double s = 0;
                while (fscanf(ptr, "%lf %lf", &tmp1, &tmp2) != EOF) {
                    tmp2 /= s;
                    if(i == 0) spectre[i++] = FunctionElement(tmp2, tmp1);
                    else spectre[i++] = FunctionElement(tmp2 + spectre[i - 1].getX(), tmp1);
                }
                break;
            }
            case FunctionSpectreOrigin::MATERIAL: {
                while (fscanf(ptr, "%lf %lf %lf", &tmp1, &tmp2, &temp) != EOF) {
                    spectre[i++] = FunctionElement(tmp1, tmp2 * k);
                }
                break;
            }
        }
	}
	fclose(ptr);
}

FunctionSpectre::FunctionSpectre(std::string filePath, FunctionSpectreOrigin origin) : FunctionSpectre(filePath, origin, 1) {
}

double FunctionSpectre::selectValue(double x) {
    double fx = 0;
    if(spectreLen > 0) {
        int i = 0;
        while(i < spectreLen && spectre[i].getX() < x) i++;
        if(i >= spectreLen) fx = spectre[spectreLen - 1].getFx();
        else if(i == 0 || spectre[i].getX() == x) fx = spectre[i].getFx();
        else fx = spectre[i - 1].getFx() + ((spectre[i].getFx() - spectre[i - 1].getFx()) / (spectre[i].getX() - spectre[i - 1].getX())) * x;
    } else {
        cout << "ERROR: empty spectre" << endl;
    }
    return fx;
}

FunctionSpectre::~FunctionSpectre() {
    delete[] spectre;
}
