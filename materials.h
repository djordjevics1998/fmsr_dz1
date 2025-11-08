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
        void multiplyXBy(double k);
        void multiplyFxBy(double k);
        ~FunctionElement();
};

/*enum FunctionSpectreOrigin {
  ENERGY,
  MATERIAL,
};*/

class FunctionSpectre {
    private:
        int spectreLen;
    protected:
        FunctionElement* spectre;
        FunctionSpectre(FunctionElement* spectre);
    public:
        FunctionSpectre();
        FunctionSpectre(int spectreLen);
        double selectValue(double x);
        int getSpectreLength();
        FunctionElement getElementAt(int i);
        /**
         * len - broj kolona u datoteci, ne ukljucujuci prvu
         * filePath - putanja do datoteke
         * invertAndNormalize - Fx i x menjaju mesta i Fx postaje kumulativna funkcija
         * k - koeficijent mnozenja Fx, kad invertAndNormalize == false
         */
        static FunctionSpectre** readMultipleFromFile(int len, std::string filePath, bool invertAndNormalize, double k);
        static FunctionSpectre* readFromFile(std::string filePath, bool invertAndNormalize, double k);
        ~FunctionSpectre();
};

#endif
