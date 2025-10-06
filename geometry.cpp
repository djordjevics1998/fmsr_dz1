#include "geometry.h"
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <set>
#include <iostream>

#define SQR(a) ((a) * (a))
#define COLL_A(a, b) (((a) - (b)) * ((a) - (b)))
#define COLL_B(c_1, v_1, c_2, v_2) (((c_1) - (c_2)) * ((v_1) - (v_2)))
using namespace std;

Vector3D::Vector3D(double x, double y, double z) {
    set(x, y, z);
}

Vector3D::Vector3D(Vector3D* v1, Vector3D* v2) : Vector3D(v2->getX() - v1->getX(), v2->getY() - v1->getY(), v2->getZ() - v1->getZ()) {
}

Vector3D::Vector3D(Vector3D* v1, Vector3D* v2, bool normalize) : Vector3D(v1, v2) {
    if (normalize) this->multiply(1 / this->len());
}

Vector3D::Vector3D(Vector3D* v) : Vector3D(v->x, v->y, v->z) {
}

double Vector3D::getX() {
    return this->x;
}

double Vector3D::getY() {
    return this->y;
}

double Vector3D::getZ() {
    return this->z;
}

void Vector3D::setX(double x) {
    this->x = x;
}

void Vector3D::setY(double y) {
    this->y = y;
}

void Vector3D::setZ(double z) {
    this->z = z;
}

void Vector3D::set(double x, double y, double z) {
    setX(x);
    setY(y);
    setZ(z);
}

void Vector3D::set(Vector3D *v) {
    set(v->x, v->y, v->z);
}

double Vector3D::len() {
    return sqrt(SQR(this->x) + SQR(this->y) + SQR(this->z));
}

void Vector3D::multiply(double m) {
    this->x *= m;
    this->y *= m;
    this->z *= m;
}

double Vector3D::scalar(Vector3D* v) {
    return this->x * v->x + this->y * v->y + this->z * v->z;
}

std::string Vector3D::toString() {
    std::ostringstream ssx, ssy, ssz;
    ssx << std::scientific << this->x;
    ssy << std::scientific << this->y;
    ssz << std::scientific << this->z;
    return "Vector3D(" + ssx.str() + ", " + ssy.str() + ", " + ssz.str() + ")";
}

Vector3D Vector3D::vector(Vector3D* v1, Vector3D* v2, bool normalize) {
    Vector3D product(v1->y * v2->z - v1->z * v2->y, v1->z * v2->x - v1->x * v2->z, v1->x * v2->y - v1->y * v2->x);
    if (normalize && (product.x != 0 || product.y != 0 || product.z != 0)) product.multiply(1 / product.len());
    return product;
}

Vector3D Vector3D::projection(Vector3D* v_orig, Vector3D* v_on) {
    Vector3D projection(v_on);
    if (projection.x != 0 || projection.y != 0 || projection.z != 0) projection.multiply(1 / projection.len());
    projection.multiply(projection.scalar(v_orig));
    return projection;
}

Vector3D Vector3D::add(Vector3D* v1, Vector3D* v2) {
    return Vector3D(v1->x + v2->x, v1->y + v2->y, v1->z + v2->z);
}

Vector3D Vector3D::sub(Vector3D* from, Vector3D* what) {
    return Vector3D(from->x - what->x, from->y - what->y, from->z - what->z);
}

double Vector3D::distance(Vector3D* v1, Vector3D* v2) {
    return sqrt(COLL_A(v1->x, v2->x) + COLL_A(v1->y, v2->y) + COLL_A(v1->z, v2->z));
}

Vector3D VECTOR3D_ZERO((double)0, 0, 0);

Particle::Particle(Vector3D *p, Vector3D *v, double e) {
    this->p = new Vector3D(p);
    this->v = new Vector3D(v);
    setE(e);
}

Particle::Particle(Vector3D *p, double e) : Particle(p, &VECTOR3D_ZERO, e) {

}

Vector3D* Particle::getP() {
    return this->p;
}

Vector3D* Particle::getV() {
    return this->v;
}

double Particle::getE() {
    return this->e;
}

void Particle::setE(double e) {
    this->e = e;
}

Particle::~Particle() {
    delete p;
    delete v;
}

Intersection::Intersection(Vector3D *p1, Vector3D *p2) {
    this->p1 = p1 == NULL ? NULL : new Vector3D(p1);
    this->p2 = p2 == NULL ? NULL : new Vector3D(p2);
}

Intersection::Intersection(Particle *p, double a, double b, double c) {
    // quadratic equation a * (x ^ 2) + b * x + c = 0
    double pom = SQR(b) - 4 * a * c;
    // ako je pom == 0 matematicki postoji resenje, medjutim, nije potrebno dalje razmatrati,
    // ne znaci nam ako ne ulazimo u objekat (ako ga samo dodirujemo)
    if((a == 0 && b == 0) || pom <= 0) {
        this->p1 = NULL;
        this->p2 = NULL;
    }
    else if(a == 0) {
        //double k = -c / b;

        // jedan koren jednacine, te nije potrebno dalje razmatrati, ne znaci nam
        // ako ne ulazimo u objekat (ako ga samo dodirujemo)
        this->p1 = NULL;
        this->p2 = NULL;
    } else {
        // korenujemo za konacno izracunavanje
        pom = sqrt(pom);
        double k1 = (-b + pom) / (2 * a);
        double k2 = (-b - pom) / (2 * a);

        // ako su oba korena manja od nule to znaci da cestica ide u suprotnom smeru od objekta,
        // te cemo postaviti da ne postoje preseci
        if(k1 < 0 && k2 < 0) {
            this->p1 = NULL;
            this->p2 = NULL;
        } else {
            Vector3D v2(p->getV());
            if(k1 >= 0) {
                Vector3D v1(p->getV());
                v1.multiply(k1);
                this->p1 = new Vector3D(Vector3D::add(p->getP(), &v1));
            }
            if(k2 >= 0) {
                Vector3D v2(p->getV());
                v2.multiply(k2);
                // ako je prvi koren jednacine negativan, to znaci da je cestica vec u telu
                if(k1 < 0) {
                    this->p1 = new Vector3D(Vector3D::add(p->getP(), &v2));
                    this->p2 = NULL;
                }
                else this->p2 = new Vector3D(Vector3D::add(p->getP(), &v2));
            } else {
                this->p2 = NULL;
            }
        }
    }
}

Vector3D* Intersection::getP1() {
    return p1;
}

Vector3D* Intersection::getP2() {
    return p2;
}

void Intersection::setP1(Vector3D* p1) {
    //if(this->p1 != NULL) {
    //    delete this->p1;
    //}
    this->p1 = p1 == NULL ? NULL : new Vector3D(p1);
}

void Intersection::setP2(Vector3D* p2) {
    //if(this->p2 != NULL) {
    //    delete this->p2;
    //}
    this->p2 = p2 == NULL ? NULL : new Vector3D(p2);
}

Intersection::~Intersection() {
    delete p1;
    delete p2;
}

PhObject::PhObject(FunctionSpectre *absU) {
    this->absU = absU;
    dose = 0;
}

FunctionSpectre* PhObject::getAbsU() {
    return absU;
}

double PhObject::getDose() {
    return dose;
}

void PhObject::setDose(double dose) {
    this->dose = dose;
}

PhObject::~PhObject() {
    delete absU;
}

Ellipsoid::Ellipsoid(Vector3D* p, double a, double b, double c, FunctionSpectre *absU) : PhObject(absU) {
    this->p = new Vector3D(p);
    this->a = a;
    this->b = b;
    this->c = c;
}

Ellipsoid::Ellipsoid(Vector3D* p, double r, FunctionSpectre *absU) : Ellipsoid(p, r, r, r, absU) {
}

Vector3D* Ellipsoid::getP() {
    return p;
}

double Ellipsoid::getA() {
    return a;
}

double Ellipsoid::getB() {
    return b;
}

double Ellipsoid::getC() {
    return c;
}

Intersection Ellipsoid::intersection(Particle* particle) {
    double ex = p->getX();
    double ey = p->getY();
    double ez = p->getZ();

    double px = particle->getP()->getX();
    double py = particle->getP()->getY();
    double pz = particle->getP()->getZ();

    double dx = px - ex;
    double dy = py - ey;
    double dz = pz - ez;

    double vx = particle->getV()->getX();
    double vy = particle->getV()->getY();
    double vz = particle->getV()->getZ();
    return Intersection(particle, SQR(vx / a) + SQR(vy / b) + SQR(vz / c),
    2 * (dx * vx / SQR(a) + dy * vy / SQR(b) + dz * vz / SQR(c)), SQR(dx / a) + SQR(dy / b) + SQR(dz / c) - 1);
}

std::string Ellipsoid::toString() {
    std::ostringstream ssd, ssx, ssy, ssz;
    ssd << std::scientific << this->dose;
    ssx << std::scientific << this->a;
    ssy << std::scientific << this->b;
    ssz << std::scientific << this->c;
    return "Ellipsoid(" + ssd.str() + ", " + this->p->toString() + ", " + ssx.str() + ", " + ssy.str() + ", " + ssz.str() + ")";
}

Ellipsoid::~Ellipsoid() {
    delete p;
}

ConeYN::ConeYN(Vector3D *p, double h, double r, FunctionSpectre *absU) : PhObject(absU) {
    this->p = new Vector3D(p);
    this->h = h;
    this->r = r;
}

Vector3D* ConeYN::getP() {
    return p;
}

double ConeYN::getH() {
    return h;
}

double ConeYN::getR() {
    return r;
}

Intersection ConeYN::intersection(Particle* particle) {
    double ex = p->getX();
    double ey = p->getY();
    double ez = p->getZ();

    double px = particle->getP()->getX();
    double py = particle->getP()->getY();
    double pz = particle->getP()->getZ();

    double dx = px - ex;
    double dy = py - ey;
    double dz = pz - ez;

    double vx = particle->getV()->getX();
    double vy = particle->getV()->getY();
    double vz = particle->getV()->getZ();

    // preseci sa omotacem beskonacne parabole, kasnije cemo proveriti bazu kupe
    auto intersection = Intersection(particle, SQR(vx) + SQR(vz) - SQR(vy * r / h),
    2 * (dx * vx + dz * vz - dy * vy * SQR(r / h)), SQR(dx) + SQR(dz) - SQR(dy * r / h) /*- 1*/);

    // ako postoji barem jedan presek sa beskonacnom parabolom
    if(intersection.getP1() != NULL || intersection.getP2() != NULL) {
        // provera da li je presek van omotaca kupe za prvo resenje
        if(intersection.getP1() != NULL) {
            auto p1 = intersection.getP1();
            if(p1->getY() >= ey || p1->getY() < ey - h) {
                intersection.setP1(intersection.getP2());
                delete p1;
                delete intersection.getP2();
                intersection.setP2(NULL);
            }
        }
        if(intersection.getP1() != NULL) {
            auto p1 = intersection.getP1();
            if(p1->getY() >= ey || p1->getY() < ey - h) {
                intersection.setP1(intersection.getP2());
                delete p1;
                delete intersection.getP2();
                intersection.setP2(NULL);
            }
        }
        // provera da li je presek van omotaca kupe za drugo resenje
        if(intersection.getP2() != NULL) {
            auto p2 = intersection.getP2();
            if(p2->getY() >= ey || p2->getY() < ey - h) {
                intersection.setP2(NULL);
                delete p2;
            }
        }
    }

    // da probamo bazu kupe kao resenje
    if(intersection.getP1() == NULL || intersection.getP2() == NULL) {
        if(vy != 0) {
            double k = (ey - h - py) / vy;
            if(k > 0) {
                double fx = px + k * vx;
                double fz = pz + k * vz;
                if(SQR(fx - ex) + SQR(fz - ez) <= SQR(r)) {
                    auto newVec = new Vector3D(fx, ey - h, fz);
                    if(intersection.getP1() == NULL) {
                        intersection.setP1(newVec);
                    } else if(intersection.getP2() == NULL) {
                        intersection.setP2(newVec);
                    }
                    delete newVec;
                }
            }
        }
    }
    return intersection;
}

std::string ConeYN::toString() {
    std::ostringstream ssd, ssh, ssr;
    ssd << std::scientific << this->dose;
    ssh << std::scientific << this->h;
    ssr << std::scientific << this->r;
    return "ConeYN(" + ssd.str() + ", " + this->p->toString() + ", " + ssh.str() + ", " + ssr.str() + ")";
}

ConeYN::~ConeYN() {
    delete p;
}

DetectorZ::DetectorZ(Vector3D *p, int xn, int yn, double cellW, double cellH) : PhObject(new FunctionSpectre()) {
    this->p = new Vector3D(p);
    this->xn = xn;
    this->yn = yn;
    this->cellW = cellW;
    this->cellH = cellH;
    this->cellCounts = new int*[this->yn];
    for(int i = 0; i < this->yn; i++) {
        this->cellCounts[i] = new int[this->xn];
    }

    for(int i = 0; i < this->yn; i++) {
        for(int j = 0; j < this->xn; j++) {
            this->cellCounts[i][j] = 0;
        }
    }
}

Vector3D* DetectorZ::getP() {
    return this->p;
}

int DetectorZ::getXN() {
    return this->xn;
}

int DetectorZ::getYN() {
    return this->yn;
}

double DetectorZ::getCellW() {
    return this->cellW;
}

double DetectorZ::getCellH() {
    return this->cellH;
}

int** DetectorZ::getCellCounts() {
    return this->cellCounts;
}

void DetectorZ::addCounter(Vector3D* pos) {
    if(pos->getZ() != p->getZ()) return;

    double ex = p->getX();
    double ey = p->getY();

    double fx = pos->getX();
    double fy = pos->getY();

    double w = xn * cellW;
    double h = yn * cellH;

    if(SQR(fx - ex) <= SQR(w / 2) && SQR(fy - ey) <= SQR(h / 2)) {
        int j = round((fx - ex + w / 2) / cellW) - 1;
        //if(j < 0) cout << j;
        //if(j < 0) {
        //    cout << fx << ", " << fy << endl;
        //}
        if(j < 0) j = 0;
        if(j >= xn) j = xn - 1;

        int i = round((fy - ey + h / 2) / cellH) - 1;
        //if(i < 0) cout << i;
        //if(i < 0) {
        //    cout << fx << ", " << fy << endl;
        //}
        if(i < 0) i = 0;
        if(i >= yn) i = yn - 1;

        this->cellCounts[i][j]++;
    }
}

Intersection DetectorZ::intersection(Particle* particle) {
    auto intersection = Intersection(NULL, NULL);

    double ex = p->getX();
    double ey = p->getY();
    double ez = p->getZ();

    double px = particle->getP()->getX();
    double py = particle->getP()->getY();
    double pz = particle->getP()->getZ();

    double vx = particle->getV()->getX();
    double vy = particle->getV()->getY();
    double vz = particle->getV()->getZ();

    double w = xn * cellW;
    double h = yn * cellH;

    if(vz != 0) {
        double k = (ez - pz) / vz;

        if(k > 0) {
            double fx = px + k * vx;
            double fy = py + k * vy;

            if(SQR(fx - ex) <= SQR(w / 2) && SQR(fy - ey) <= SQR(h / 2)) {
                auto newVec = new Vector3D(fx, fy, ez);
                //cout << fx << ", " << fy << ", " << ez << endl;
                intersection.setP1(newVec);
                delete newVec;
            }
        }
    }

    return intersection;
}

std::string DetectorZ::toString() {
    std::ostringstream ssCellH, ssCellW;
    ssCellW << std::scientific << this->cellW;
    ssCellH << std::scientific << this->cellH;
    return "DetectorZ(" + this->p->toString() + ", " + std::to_string(xn) + ", " + std::to_string(yn) + ", " + ssCellW.str() + ", " + ssCellH.str() + ")";
}

DetectorZ::~DetectorZ() {
    delete p;
    for(int i = 0; i < yn; i++) {
        delete[] cellCounts[i];
    }
    delete[] cellCounts;
}

Event::Event(Particle *particle, PhObject *object, Intersection *intersection) {
    this->particle = particle;
    this->object = object;
    this->intersection = intersection;
}

Particle* Event::getParticle() {
    return this->particle;
}

PhObject* Event::getObject() {
    return this->object;
}

Intersection* Event::getIntersection() {
    return this->intersection;
}

double Event::getMinDistance() {
    double distance = -1;
    if(intersection->getP1() != NULL || intersection->getP2() != NULL) {
        if(intersection->getP1() != NULL) {
            double cand = Vector3D::distance(particle->getP(), intersection->getP1());
            //cout << cand << "|" << intersection->getP1()->toString() << "|";
            if(cand < distance || distance == -1) distance = cand;
        }

        if(intersection->getP2() != NULL) {
            double cand = Vector3D::distance(particle->getP(), intersection->getP2());
            //cout << cand << "|" << intersection->getP2()->toString() << "|";
            if(cand < distance || distance == -1) distance = cand;
        }
    }
    //cout << distance << endl;
    //cout << Vector3D::distance(intersection->getP1(), intersection->getP2()) << endl;
    return distance;
}

bool Event::compare(Event* o1, Event* o2) {
    return o1->getMinDistance() <= o2->getMinDistance();
}

Vector3D *Event::getCloserIntersection() {
    if(intersection->getP1() == NULL) return intersection->getP2();
    if(intersection->getP2() == NULL) return intersection->getP1();
    return Vector3D::distance(particle->getP(), intersection->getP1()) <=
        Vector3D::distance(particle->getP(), intersection->getP2()) ? intersection->getP1() : intersection->getP2();
}

Vector3D *Event::getFartherIntersection() {
    if(intersection->getP1() == NULL) return intersection->getP2();
    if(intersection->getP2() == NULL) return intersection->getP1();
    return Vector3D::distance(particle->getP(), intersection->getP1()) >=
        Vector3D::distance(particle->getP(), intersection->getP2()) ? intersection->getP1() : intersection->getP2();
}

Event::~Event() {
    delete intersection;
}
