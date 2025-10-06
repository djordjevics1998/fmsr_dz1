#ifndef H_MATERIALS
#define H_MATERIALS

#include <string>

class FunctionElement {
    protected:
        double x, fx;
    public:
        FunctionElement();
        FunctionElement(double x, double fx);
        double getX();
        double getFx();
        ~FunctionElement();
};

enum FunctionSpectreOrigin {
  ENERGY,
  MATERIAL,
};

class FunctionSpectre {
    private:
        int spectreLen;
    protected:
        FunctionElement* spectre;
    public:
        FunctionSpectre();
        FunctionSpectre(std::string filePath, FunctionSpectreOrigin origin, double k);
        FunctionSpectre(std::string filePath, FunctionSpectreOrigin origin);
        double selectValue(double x);
        ~FunctionSpectre();
};

#endif
