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

int main(int argc, char *argv[]) {
	mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
	uniform_real_distribution<> distR(0, 1);
	double limitE = 10; // keV, ispod koje se foton vise ne razmatra

    cout << "Unesite stepen broja istorija (10 ^ X):";
    int stepen;
    cin >> stepen;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    long long N = 1;
    int i, k;
    for(i = 0; i < stepen; i++) N *= 10;

    FunctionSpectre* energySpec = FunctionSpectre::readFromFile("120kVp 1mm.txt", true, 1);

    auto particlePos0 = new Vector3D((double)0, 0, -1000);
    // "Water, Liquid.txt"
    auto phantom = new Ellipsoid(new Vector3D((double) 0, 0, 0), 300, 300, 60, FunctionSpectre::readMultipleFromFile(2, "voda1.txt", false, 1.0), 2);
    auto detector = new DetectorZ(new Vector3D((double)0, 0, 70), 380, 380, 1, 1);
    PhObject* objects[] = {
        // "Bone, Cortical.txt"
        new Ellipsoid(new Vector3D((double) 0, -100, 0), 40, 40, 40, FunctionSpectre::readMultipleFromFile(2, "kost1.txt", false, 1.92), 2),
        // "Blood, Whole.txt"
        new ConeYN(new Vector3D((double)0, 150, 0), 120, 30, FunctionSpectre::readMultipleFromFile(2, "krvjod1.txt", false, 1.06), 2),
    };
    int objectsLen = sizeof(objects) / sizeof(PhObject*);

    auto kLinije = FunctionSpectre::readMultipleFromFile(1, "JodKlinije.txt", false, 1000);
    double sumProb = 0;
    for(i = 0; i < kLinije[0]->getSpectreLength(); i++) {
        kLinije[0]->getElementAt(i).multiplyXBy(1 / 1000.0);
        sumProb += kLinije[0]->getElementAt(i).getX();
    }

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
        auto particle = new Particle(particlePos0, new Vector3D(sin(theta) * cos(phi),  sin(theta) * sin(phi), cos(theta)), energySpec->selectValue(distR(rng)));
        //cout << theta << " ";
        //cout << particle->getV()->getZ() << endl;
        auto phin = phantom->intersection(particle);

        if(phin.getP1() != NULL || phin.getP2() != NULL) {
            N--;

            bool shouldContinue = false;
            bool shouldBreak = false;
            bool isKLinija = false;
            // petlja kretanja fotona
            while(particle->getE() >= limitE && !shouldBreak) {
                shouldBreak = false;

                // odredjivanje u okviru kojeg je objekta foton trenutno
                PhObject* insideObj = NULL;
                /*if(isScattered) {
                    cout << particle->getV()->getZ() << endl;
                }*/

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
                    //if(sectsLen > 3) cout << origSectsLen << " " << sectsLen << endl;
                } else if(phantomIntersection.getP1() != NULL && phantomIntersection.getP2() != NULL) {
                    auto phantomEvent = new Event(particle, phantom, new Intersection(phantomIntersection.getP1(), phantomIntersection.getP2()));
                    events[sectsLen++] = new Event(particle, phantom, new Intersection(phantomEvent->getCloserIntersection(), phantomEvent->getFartherIntersection()));
                    delete phantomEvent;
                }

                double domet = -log(distR(rng));
                shouldContinue = false;
                for(i = 0; i < sectsLen; i++) {
                    auto fsu = events[i]->getObject()->getFsU();
                    auto e = events[i]->getParticle()->getE();
                    // delimo sa 10 da bi dobili cm
                    double mult = Vector3D::distance(events[i]->getIntersection()->getP1(), events[i]->getIntersection()->getP2()) / 10;
                    //double dtabs = fsu[1]->selectValue(e) * mult;
                    //double dtsca = fsu[0]->selectValue(e) * mult;

                    double uabs = fsu[1]->selectValue(e);
                    double uras = fsu[0]->selectValue(e);
                    double duk = (uabs + uras) * mult;
                    if(domet <= duk) {
                        if(distR(rng) <= uabs / (uabs + uras)) {
                            // apsorpcija
                            shouldBreak = true;

                            // objects[1] je sa k jodom
                            if(!shouldContinue && !isKLinija && events[i]->getObject() == objects[1] && distR(rng) <= 0.33) {
                                auto particle = events[i]->getParticle();
                                particle->getV()->set(distR(rng), distR(rng), distR(rng));
                                if(particle->getV()->len() != 0) particle->getV()->multiply(1 / particle->getV()->len());

                                auto oldE = particle->getE();
                                auto prob = distR(rng) * sumProb;
                                double cumProb = 0;
                                int ib = 0;
                                while(ib < kLinije[0]->getSpectreLength() && prob > cumProb + kLinije[0]->getElementAt(ib).getX()) {
                                    cumProb += kLinije[0]->getElementAt(ib).getX();
                                    ib++;
                                }
                                if(ib >= kLinije[0]->getSpectreLength()) {
                                    ib = kLinije[0]->getSpectreLength() - 1;
                                }

                                if(kLinije[0]->getElementAt(ib).getFx() > 0 && oldE > kLinije[0]->getElementAt(ib).getFx()) {
                                    shouldContinue = true;
                                    isKLinija = true;
                                    shouldBreak = false;

                                    particle->setE(kLinije[0]->getElementAt(ib).getFx());
                                    events[i]->getObject()->setDose(events[i]->getObject()->getDose() + (oldE - events[i]->getParticle()->getE()));
                                }
                            }

                            if(!shouldContinue) events[i]->getObject()->setDose(events[i]->getObject()->getDose() + events[i]->getParticle()->getE());
                        } else {
                            shouldBreak = false;
                            shouldContinue = true;
                            if(true) {//events[i]->getParticle()->getV()->len() != 0) {
                                // mnozimo sa 10 jer smo malopre delili sa 10
                                double k = 10 * (domet / (uabs + uras)); // / events[i]->getParticle()->getV()->len();
                                //cout << k << " ";
                                events[i]->getParticle()->getP()->set(events[i]->getCloserIntersection());
                                events[i]->getParticle()->getP()->set(events[i]->getParticle()->getP()->getX() + events[i]->getParticle()->getV()->getX() * k,
                                    events[i]->getParticle()->getP()->getY() + events[i]->getParticle()->getV()->getY() * k,
                                    events[i]->getParticle()->getP()->getZ() + events[i]->getParticle()->getV()->getZ() * k
                                );
                                //cout << events[i]->getParticle()->getP()->getX() << " " << events[i]->getParticle()->getP()->getY() << " " << events[i]->getParticle()->getP()->getZ() << endl;
                            }
                            double olde = events[i]->getParticle()->getE();
                            events[i]->getParticle()->scatter(rng);
                            events[i]->getObject()->setDose(events[i]->getObject()->getDose() + olde - events[i]->getParticle()->getE());

                            // bice apsorbovano
                            if(particle->getE() < limitE)
                                events[i]->getObject()->setDose(events[i]->getObject()->getDose() + particle->getE());
                        }
                        break;
                    }
                    domet -= duk;

                    /*if(domet <= dtabs) {
                        shouldBreak = true;
                        events[i]->getObject()->setDose(events[i]->getObject()->getDose() + events[i]->getParticle()->getE());
                        break;
                    }
                    domet -= dtabs;
                    if(domet <= dtsca) {
                        shouldBreak = false;
                        isScattered = true;
                        if(events[i]->getParticle()->getV()->len() != 0) {
                            double k = (domet / fsu[0]->selectValue(e)) / events[i]->getParticle()->getV()->len();
                            events[i]->getParticle()->getP()->set(events[i]->getParticle()->getP()->getX() + events[i]->getParticle()->getV()->getX() * k,
                                events[i]->getParticle()->getP()->getY() + events[i]->getParticle()->getV()->getY() * k,
                                events[i]->getParticle()->getP()->getZ() + events[i]->getParticle()->getV()->getZ() * k
                            );
                        }
                        double olde = events[i]->getParticle()->getE();
                        events[i]->getParticle()->scatter(rng);
                        events[i]->getObject()->setDose(events[i]->getObject()->getDose() + olde - events[i]->getParticle()->getE());
                        break;
                    }
                    domet -= dtsca;*/
                }
                //shouldBreak = sectsLen < 2;

                if(!shouldBreak && !shouldContinue) {
                    auto intersection = detector->intersection(particle);
                    if(intersection.getP1() != NULL) {
                        auto pos = intersection.getP1();
                        detector->addCounter(pos);
                    }
                    shouldBreak = true;
                }

                for(i = 0; i < sectsLen; i++) delete events[i];
            }
        } else {
            cout << "PROMASIO" << endl;
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
    long long uhv = 0;
	for(i = 0; i < detector->getYN(); i++) {
		for(int k = 0; k < detector->getXN(); k++) {
            uhv += detector->getCellCounts()[i][k];
            fprintf(ptr, "%d ", (int) (4095 * (detector->getCellCounts()[i][k] - mind) / (maxd - mind)));
		}
		fprintf(ptr, "\n");
	}
	cout << "Broj uhvacenih fotona: " << uhv << endl;
	fclose(ptr);
	printf("Simulacija je zavrsena. Slika je sacuvana kao cpp_slika.txt");

	delete phantom;
	delete detector;
    for(i = 0; i < objectsLen; i++) {
        delete objects[i];
    }
    delete energySpec;
    delete kLinije[0];
    delete kLinije;

    //µ
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << endl << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    return 0;
}
