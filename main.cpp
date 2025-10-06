#include "materials.h"
#include "geometry.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <random>
#include <chrono>
#include <algorithm>
#include <iostream>

#define SQR(a) ((a) * (a))

using namespace std;

struct Element {
	double ro;
	char name[30];
} WATER_LIQUID = {1.0f, "Water, Liquid"}, BLOOD_WHOLE = {1.06, "Blood, Whole"}, BONE_CORTICAL = {1.92, "Bone, Cortical"};

int main(int argc, char *argv[]) {
	mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
	uniform_real_distribution<> distR(0, 1);
	double limitE = 10; // keV, ispod koje se foton vise ne razmatra

    cout << "Unesite stepen broja istorija (10 ^ X):";
    int stepen;
    cin >> stepen;
    long long N = 1;
    int i, k;
    for(i = 0; i < stepen; i++) N *= 10;

    FunctionSpectre energySpec("120kVp 1mm.txt", FunctionSpectreOrigin::ENERGY);

    auto particlePos0 = new Vector3D((double)0, 0, -1000);
    auto phantom = new Ellipsoid(new Vector3D((double) 0, 0, 0), 300, 300, 60, new FunctionSpectre("Water, Liquid.txt", FunctionSpectreOrigin::MATERIAL, 1.0));
    auto detector = new DetectorZ(new Vector3D((double)0, 0, 70), 380, 380, 1, 1);
    PhObject* objects[] = {
        new Ellipsoid(new Vector3D((double) 0, -100, 0), 40, 40, 40, new FunctionSpectre("Bone, Cortical.txt", FunctionSpectreOrigin::MATERIAL, 1.92)),
        new ConeYN(new Vector3D((double)0, 150, 0), 120, 30, new FunctionSpectre("Blood, Whole.txt", FunctionSpectreOrigin::MATERIAL, 1.06)),
    };
    int objectsLen = sizeof(objects) / sizeof(PhObject*);
    // maksimalan broj trajektorija: svako telo u fantomu,
    // izmedju svaka dva tela, pre prvog tela u fantomu i
    // posle poslednjeg tela u fantomu
    int maxSectsLen = objectsLen * 2 + 1;
    Event* events[maxSectsLen];
    auto maxR = max(phantom->getA(), phantom->getB());
    double thetaGr = atan(maxR / abs(particlePos0->getZ() - phantom->getP()->getZ()));
    while(N > 0) {
        //auto theta = thetaGr * distR(rng);
        auto theta = acos(1 - (1 - cos(thetaGr)) * distR(rng));
        auto phi = 2 * M_PI * distR(rng);

        //auto particle = new Particle(particlePos0, new Vector3D(sin(theta) * cos(phi),  sin(theta) * sin(phi), cos(theta)), e);
        auto particle = new Particle(particlePos0, new Vector3D(tan(theta) * cos(phi),  tan(theta) * sin(phi), 1 - tan(theta)), energySpec.selectValue(distR(rng)));

        auto phin = phantom->intersection(particle);

        if(phin.getP1() != NULL || phin.getP2() != NULL) {
            N--;

            bool shouldBreak = false;
            // petlja kretanja fotona
            while(particle->getE() >= limitE && !shouldBreak) {
                shouldBreak = false;

                // odredjivanje u okviru kojeg je objekta foton trenutno
                PhObject* insideObj = NULL;

                int sectsLen = 0;
                for(i = 0; i < objectsLen; i++) {
                    auto intersection = objects[i]->intersection(particle);
                    // ne interesuje nas objekat kroz koji ne prolazimo
                    if(intersection.getP1() == NULL && intersection.getP2() == NULL) {
                        //delete intersection;
                        continue;
                    }
                    if(intersection.getP1() == NULL) {
                        insideObj = objects[i];
                        intersection.setP1(new Vector3D(particle->getP()));
                    } else if(intersection.getP2() == NULL) {
                        insideObj = objects[i];
                        intersection.setP2(intersection.getP1());
                        intersection.setP1(new Vector3D(particle->getP()));
                    }
                    //cout << intersection.getP1()->toString() << ", " << intersection.getP2()->toString() << endl;
                    events[sectsLen++] = new Event(particle, objects[i], new Intersection(intersection.getP1(), intersection.getP2()));
                }
                auto phantomIntersection = phantom->intersection(particle);
                if(insideObj == NULL) {
                    if(phantomIntersection.getP1() != NULL || phantomIntersection.getP2() != NULL) {
                        if(phantomIntersection.getP1() == NULL) {
                            insideObj = phantom;
                            phantomIntersection.setP1(new Vector3D(particle->getP()));
                        } else if(phantomIntersection.getP2() == NULL) {
                            insideObj = phantom;
                            phantomIntersection.setP2(phantomIntersection.getP1());
                            phantomIntersection.setP1(new Vector3D(particle->getP()));
                        }
                    }
                }
                if(sectsLen > 0) {
                    sort(events, events + sectsLen, Event::compare);
                    /*if(sectsLen > 1) {
                        cout << events[0]->getObject()->toString() << ", " << events[1]->getObject()->toString() << endl;
                        cout << events[0]->getMinDistance() << ", " << events[1]->getMinDistance() << endl;
                        cout << events[0]->getObject()->toString() << endl;
                    }*/
                    int origSectsLen = sectsLen;
                    auto phantomEvent = new Event(particle, phantom, new Intersection(phantomIntersection.getP1(), phantomIntersection.getP2()));
                    if(insideObj == NULL || insideObj == phantom) {
                        if(insideObj == phantom) {
                            events[sectsLen++] = new Event(particle, phantom, new Intersection(particle->getP(), events[0]->getCloserIntersection()));
                        } else if(insideObj == NULL) {
                            events[sectsLen++] = new Event(particle, phantom, new Intersection(phantomEvent->getCloserIntersection(), events[0]->getCloserIntersection()));
                        }
                    }
                    for(i = 1; i < origSectsLen; i++) {
                        events[sectsLen++] = new Event(particle, phantom, new Intersection(events[i - 1]->getFartherIntersection(), events[i]->getCloserIntersection()));
                    }
                    events[sectsLen++] = new Event(particle, phantom, new Intersection(events[origSectsLen - 1]->getFartherIntersection(), phantomEvent->getFartherIntersection()));
                    delete phantomEvent;
                    sort(events, events + sectsLen, Event::compare);
                } else if(phantomIntersection.getP1() != NULL && phantomIntersection.getP2() != NULL) {
                    auto phantomEvent = new Event(particle, phantom, new Intersection(phantomIntersection.getP1(), phantomIntersection.getP2()));
                    events[sectsLen++] = new Event(particle, phantom, new Intersection(phantomEvent->getCloserIntersection(), phantomEvent->getFartherIntersection()));
                    delete phantomEvent;
                }

                double domet = -log(distR(rng));
                for(i = 0; i < sectsLen; i++) {
                    double dt = events[i]->getObject()->getAbsU()->selectValue(events[i]->getParticle()->getE())
                     * Vector3D::distance(events[i]->getIntersection()->getP1(), events[i]->getIntersection()->getP2());
                    if(domet <= dt) {
                        shouldBreak = true;
                        break;
                    }
                    domet -= dt;
                }
                //shouldBreak = sectsLen < 2;

                if(!shouldBreak) {
                    auto intersection = detector->intersection(particle);
                    if(intersection.getP1() != NULL) {
                        auto pos = intersection.getP1();
                        detector->addCounter(pos);
                    }
                    shouldBreak = true;
                }

                for(i = 0; i < sectsLen; i++) delete events[i];
                shouldBreak = true;
            }
        } else {
        }

        delete particle;
    }

    auto img = detector->getCellCounts();
    double mind = img[0][0], maxd = img[0][0];
	for(i = 0; i < detector->getYN(); i++)
		for(k = 0; k < detector->getXN(); k++) {
			if(img[i][k] > maxd) maxd = img[i][k];
			if(img[i][k] < mind) mind = img[i][k];
		}

    auto ptr = fopen("cpp_slika.txt", "w");
	for(i = 0; i < detector->getYN(); i++) {
		for(int k = 0; k < detector->getXN(); k++) fprintf(ptr, "%d ", (int) (4095 * (detector->getCellCounts()[i][k] - mind) / (maxd - mind)));
		fprintf(ptr, "\n");
	}
	fclose(ptr);
	printf("Simulacija je zavrsena. Slika je sacuvana kao cpp_slika.txt");

	delete phantom;
	delete detector;
    for(i = 0; i < objectsLen; i++) {
        delete objects[i];
    }

    return 0;
}
