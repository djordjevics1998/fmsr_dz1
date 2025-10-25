#ifndef H_GEOMETRY
#define H_GEOMETRY

#include <string>
#include <random>
#include "materials.h"

class Vector3D {
protected:
	double x, y, z;

public:
	Vector3D(double x, double y, double z);
	Vector3D(Vector3D* v1, Vector3D* v2);
	Vector3D(Vector3D* v1, Vector3D* v2, bool normalize);
	Vector3D(Vector3D* v);

	double getX();
	double getY();
	double getZ();
	void setX(double x);
	void setY(double y);
	void setZ(double z);
	void set(double x, double y, double z);
	void set(Vector3D* v);
	double len();
	void multiply(double m);
	double scalar(Vector3D* v);
	std::string toString();

	static Vector3D vector(Vector3D* v1, Vector3D* v2, bool normalize);
	static Vector3D projection(Vector3D* v_orig, Vector3D* v_on);
	static Vector3D add(Vector3D* v1, Vector3D* v2);
	static Vector3D sub(Vector3D* from, Vector3D* what);
	static double distance(Vector3D* v1, Vector3D* v2);
};

class Particle {
	protected:
		Vector3D *p, *v;
		double e;

	public:
		Particle(Vector3D *p, Vector3D *v, double e);
		Particle(Vector3D *p, double e);
		Vector3D* getP();
		Vector3D* getV();
		double getE();
		void setE(double e);
		void scatter(std::mt19937 rng);
		~Particle();
};

class Intersection {
	protected:
		Vector3D *p1, *p2;
	public:
		Intersection(Vector3D *p1, Vector3D *p2);
		/**
		 * @return vrednost je uvek objekat Intersection cija svojstva ako su razlicite od NULL oznacavaju
		 * poziciju/e preseka koje ce cestica sigurno pogoditi. Ako je cestica npr. vec u telu, onda ce samo
		 * Intersection::getP1() biti razlicito od NULL i oznaciti tacku gde cestica izlazi iz tela
		 */
		Intersection(Particle *p, double a, double b, double c);
		Vector3D* getP1();
		Vector3D* getP2();
		void setP1(Vector3D* p1);
		void setP2(Vector3D* p2);
		~Intersection();
};

class PhObject {
	protected:
        FunctionSpectre **fsu;
        int fsuLen;
		double dose;
	public:
		PhObject(FunctionSpectre **fsu, int fsuLen);
		virtual Intersection intersection(Particle* particle) = 0;
		FunctionSpectre** getFsU();
		double getDose();
		void setDose(double dose);
		virtual std::string toString() = 0;
		virtual ~PhObject();
};

class Ellipsoid : public PhObject {
	protected:
		Vector3D* p;
		double a, b, c;
	public:
		Ellipsoid(Vector3D* p, double a, double b, double c, FunctionSpectre **absU, int fsuLen);
		/**
		 * Sphere
		 */
		Ellipsoid(Vector3D* p, double r, FunctionSpectre **absU, int fsuLen);
		Vector3D* getP();
		double getA();
		double getB();
		double getC();
        Intersection intersection(Particle *particle) override;
		std::string toString() override;
        ~Ellipsoid();
};

/**
 * Kupa sa osom paralelnom Y osi, okrenuta nadole (ka negativnom kraju Y ose)
 */
class ConeYN : public PhObject {
	protected:
		Vector3D *p;
		double h, r;
	public:
		ConeYN(Vector3D *p, double h, double r, FunctionSpectre **absU, int fsuLen);
		Vector3D* getP();
		double getH();
		double getR();
        Intersection intersection(Particle *particle) override;
		std::string toString() override;
		~ConeYN();
};

class DetectorZ : public PhObject {
	protected:
		Vector3D *p;
		int xn, yn;
		double cellW, cellH;
		int** cellCounts;
	public:
		DetectorZ(Vector3D *p, int xn, int yn, double cellW, double cellH);
		Vector3D* getP();
		int getXN();
		int getYN();
		double getCellW();
		double getCellH();
		int** getCellCounts();
		void addCounter(Vector3D *pos);
        Intersection intersection(Particle *particle) override;
		std::string toString() override;
		~DetectorZ();
};

class Event {
	protected:
		Particle* particle;
		PhObject* object;
		Intersection* intersection;
	public:
		Event(Particle* particle, PhObject* object, Intersection* intersection);
		Particle* getParticle();
		PhObject* getObject();
		Intersection* getIntersection();
		double getMinDistance();
		Vector3D* getCloserIntersection();
		Vector3D* getFartherIntersection();
		static bool compare(Event* e1, Event* e2);
		~Event();
};

#endif
