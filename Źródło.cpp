#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fstream>

using namespace std;

#define schematCalkowania 3
#define size schematCalkowania*schematCalkowania

class Element
{
public:
	int ID[4];
	double H[4][4]; //H lokalne
	double C[4][4]; //C lokalne
	//double Hbc[4][4];
	double P[4]; //P local 
};Element* Elem;

class Node
{
public:
	//wspó³rzedne wêz³ów
	double x;
	double y;
	//dla C
	double t0 = 100; //temperatura pocz¹tkowa 
	//dla Hbc
	int Bc = 0; //warunek brzegowy 
}; Node* Nd;

class Data
{
public:
	/*H- wysokosc siatki,
	W- szerokoœæ siatki,
	nH- liczba wezlow po wysokosci
	nW- liczba wezlow po szerokosci
	nE- liczba elementow
	nN- liczba wezlow
	k, ro, cp -parametry materia³owe
	*/
	double t0, timeS, deltaTau, Temp_otoczenia, alfa, H, W, nH, nW, c, k, ro;
	double nE, nN, dX, dY;
	//Danie do macierzy H:
	//double H = 0.1;
//	double W = 0.1;
	//int nH = 4;
//	int nW = 4;
	//double k = 25.0;  //wpó³czynik przewodzenia ciep³a [W/(m*C)]
	//Dane dla macierzy C:
	//double t0 = 100;  //temperatura pocz¹tkowa
	//double c = 700; //cieplo wlasciwe [J/(kg*C)]
	//double ro = 7800; //gestosc [kg/m3]
	//Dane dla macierzy Hbc:
	//double alfa = 300;  //[W/m2K]
	//do P 
	//double Temp_otoczenia = 1200; //temp. otoczenia
	//do H zastepczej
	//double deltaTau = 50; //[s]
	//double timeS = 500; //[s]

	Data() {
		fstream plik;
		plik.open("test1.txt", ios::in);
		plik >> t0; //100- temperatura pocz¹tkowa [C]
		plik >> timeS; //500s- ca³kowity czas rozwi¹zania [Hz]
		plik >> deltaTau; //50s- krok czasowy (czas jednej iteracji)  [Hz]
		plik >> Temp_otoczenia; //1200
		plik >> alfa; //300- (W/m2K) konwekcyna wymiana ciep³a [Hbc] [P]
		plik >> H; //0.1 [m]
		plik >> W; //0.1 [m]
		plik >> nH; //4 
		plik >> nW; //4
		plik >> c; //700 [J/kgC]- ciep³o w³aœciwe  [C]
		plik >> k; //25 [W/mC]- wspó³æzynnik przewodzenia ciep³a [H]
		plik >> ro; //7800 [kg/m3]- gêstoœæ [C]
		plik.close(); 
		 nE = (nH - 1) * (nW - 1); //3*3
		 nN = nH * nW;
		 dX = W / (nW - 1);
		 dY = H / (nH - 1);
	}

}; Data* Dt = new Data();

void creatingMesh() {
	Elem = new Element[Dt->nE];
	Nd = new Node[Dt->nN];

	//WSPÓ£RZÊDNE WÊZ£ÓW
	int iNode = 0; //iterator wêz³owy

	//cout << "Wspolrzedne wezlow: " << endl;
	for (int i = 0; i < Dt->nW; i++)
	{
		double y0 = 0;
		for (int j = 0; j < Dt->nH; j++)
		{
			Nd[iNode].x = i * Dt->dX;
			Nd[iNode].y = y0 + j * Dt->dY;

			if (Nd[iNode].x == 0 || Nd[iNode].x == Dt->W || Nd[iNode].y == 0 || Nd[iNode].y == Dt->H) {
				Nd[iNode].Bc = 1;
			}

			//	cout << "Wezel nr. " << iNode << " (x,y)=(" << Nd[iNode].x << "," << Nd[iNode].y << ")" << endl;
			iNode++;
		}
	}



	//	cout << endl << "Id wezlow: " << endl;
		//int k = 0;
	for (int i = 0; i < Dt->nE; i++)
	{
		//if ((k % Dt->nH) == 0) k++;
		Elem[i].ID[0] = i + i / (Dt->nH - 1);         // k;
		Elem[i].ID[1] = Elem[i].ID[0] + Dt->nH;
		Elem[i].ID[2] = Elem[i].ID[1] + 1;
		Elem[i].ID[3] = Elem[i].ID[0] + 1;
		//k++;
		//cout << "Dla elementu: " << i << " ID[" << Elem[i].ID[0] << ", " << Elem[i].ID[1] << ", " << Elem[i].ID[2] << " ," << Elem[i].ID[3] << "]" << endl;
	}
}

class Derivatives
{
public: 
	//eta(os Y) ksi(os X)
	double N[size][4];
	double dNdEta[size][4];
	double dNdKsi[size][4];	
	
	/*
	//2 pkt schemat ca³kowania
	double eta[size] = { -1/sqrt(3), -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3) };
	double ksi[size] = { -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3) };
	double w1[size] = { 1, 1, 1, 1 };
	double w2[size] = { 1,1,1,1 };
	double ksiHbc[schematCalkowania * 4] = { -1 / sqrt(3),  1 / sqrt(3), 1, 1, 1 / sqrt(3),  -1 / sqrt(3), -1, -1 };
	double etaHbc[schematCalkowania * 4] = { -1, -1, -1 / sqrt(3),  1 / sqrt(3), 1 , 1, 1 / sqrt(3),  -1 / sqrt(3) };
	double wagaHbc[schematCalkowania] = { 1, 1 };
	double Nbc[schematCalkowania * 4][4];
	*/
	
	
	//3 pkt schemat calkowania
	double ksi[size] = { -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0) };
	double w2[size] = { 5.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 5.0/9.0 };
	double eta[size] = { -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, 0, 0, sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), sqrt(3.0 / 5.0) };
	double w1[size] = { 5.0/9.0, 5.0/9.0, 5.0/9.0, 8.0/9.0, 8.0/9.0, 8.0/9.0, 5.0/9.0, 5.0/9.0, 5.0/9.0 };

	double ksiHbc[schematCalkowania * 4] = { -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), 1, 1, 1, sqrt(3.0 / 5.0), 0, -sqrt(3.0 / 5.0), -1, -1, -1 };
	double etaHbc[schematCalkowania * 4] = { -1, -1, -1, -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0), 1, 1, 1, sqrt(3.0 / 5.0), 0, -sqrt(3.0 / 5.0) };
	double wagaHbc[schematCalkowania] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	double Nbc[schematCalkowania * 4][4];
	
	
	/*
	//4 pkt schemat calkowania
	double ksi[size] = { -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) };
	double w2[size] = { 0.347855, 0.347855, 0.347855, 0.347855, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.652145, 0.6521450, 0.652145,0.347855, 0.347855, 0.347855, 0.347855 };
	double w1[size] = { 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855, 0.347855, 0.652145, 0.652145, 0.347855 };
	double eta[size] = { -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)),sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)),sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)),sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) };
	
	double ksiHbc[schematCalkowania*4]= { -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), 1, 1, 1, 1, sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -1, -1, -1, -1 };
	double etaHbc[schematCalkowania*4]= { -1, -1, -1, -1, -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), 1, 1, 1, 1, sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) };
	double wagaHbc[schematCalkowania] = { (18.0 - sqrt(30)) / 36.0, (18.0 + sqrt(30)) / 36.0, (18.0 + sqrt(30)) / 36.0, (18.0 - sqrt(30)) / 36.0 };
	double Nbc[schematCalkowania * 4][4];
	*/
	Derivatives() {
			//dla macierzy Hbc
			for (int i = 0; i < schematCalkowania * 4; i++)
			{
				//if (i % 3 == 0)// cout << endl;
				//[pktCA£KOWAIA][POWIERZCHNIA]
				Nbc[i][0] = 0.25 * ((1 - ksiHbc[i]) * (1 - etaHbc[i]));
				Nbc[i][1] = 0.25 * ((1 + ksiHbc[i]) * (1 - etaHbc[i]));
				Nbc[i][2] = 0.25 * ((1 + ksiHbc[i]) * (1 + etaHbc[i]));
				Nbc[i][3] = 0.25 * ((1 - ksiHbc[i]) * (1 + etaHbc[i]));
				//	cout << Nbc[i][0] << ",   " << Nbc[i][1] << ",   " << Nbc[i][2] << ",   " << Nbc[i][3] << endl;

			}

			//dla macierzy C
			for (int i = 0; i < size; i++) {
				N[i][0] = (0.25) * (1 - ksi[i]) * (1 - eta[i]);
				N[i][1] = (0.25) * (1 + ksi[i]) * (1 - eta[i]);
				N[i][2] = (0.25) * (1 + ksi[i]) * (1 + eta[i]);
				N[i][3] = (0.25) * (1 - ksi[i]) * (1 + eta[i]);
			//cout << "{" << N[i][0] << ", " << N[i][1] << ", " << N[i][2] << ", " << N[i][3] << "} " << endl;
			
			}
			
			//dla macierzy H
			for (int i = 0; i < size; i++) {
				dNdEta[i][0] = (-0.25) * (1 - ksi[i]);
				dNdEta[i][1] = (-0.25) * (1 + ksi[i]);
				dNdEta[i][2] = (0.25) * (1 + ksi[i]);
				dNdEta[i][3] = (0.25) * (1 - ksi[i]);
				//cout << "{" << dNdEta[i][0] << ", " << dNdEta[i][1] << ", " << dNdEta[i][2] << ", " << dNdEta[i][3] << "} " << endl;
				
			}
			for (int i = 0; i < size; i++) {
				dNdKsi[i][0] = (-0.25) * (1 - eta[i]);
				dNdKsi[i][1] = (0.25) * (1 - eta[i]);
				dNdKsi[i][2] = (0.25) * (1 + eta[i]);
				dNdKsi[i][3] = (-0.25) * (1 + eta[i]);
				//	cout << "{" << dNdKsi[i][0] << ", " << dNdKsi[i][1] << ", " << dNdKsi[i][2] << ", " << dNdKsi[i][3] << "} " << endl;
					
			}
		
	}
}; Derivatives* el4;

double createJacobian(double** jakobian, double* x, double* y, int pkt, double** odwrJ)
{
	//INTERPOLACJA WSPÓ£RZÊDNEJ X 
	//   dx/dksi
	jakobian[0][0] = el4->dNdKsi[pkt][0] * x[0] + el4->dNdKsi[pkt][1] * x[1] + el4->dNdKsi[pkt][2] * x[2] + el4->dNdKsi[pkt][3] * x[3];
	//dx/deta
	jakobian[1][0] = el4->dNdEta[pkt][0] * x[0] + el4->dNdEta[pkt][1] * x[1] + el4->dNdEta[pkt][2] * x[2] + el4->dNdEta[pkt][3] * x[3];
	//INTERPOLACJA WSPÓ£RZÊDNEJ Y
	//   dy/dksi
	jakobian[0][1] = el4->dNdKsi[pkt][0] * y[0] + el4->dNdKsi[pkt][1] * y[1] + el4->dNdKsi[pkt][2] * y[2] + el4->dNdKsi[pkt][3] * y[3];
	//   dy/deta
	jakobian[1][1] = el4->dNdEta[pkt][0] * y[0] + el4->dNdEta[pkt][1] * y[1] + el4->dNdEta[pkt][2] * y[2] + el4->dNdEta[pkt][3] * y[3];

	//cout << endl << "Jakobian dla " << pkt << " punktu ca³kowania" << endl << endl;
	//cout << jakobian[0][0] << "  " << jakobian[0][1] << endl << jakobian[1][0] << "  " << jakobian[1][1] << endl << endl;

	//odwrotnosc jakobianow 
	odwrJ[0][0] = jakobian[1][1];  odwrJ[0][1] = -jakobian[0][1];
	odwrJ[1][0] = -jakobian[1][0]; odwrJ[1][1] = jakobian[0][0];

	//cout << endl << "Odwrtnosc jakobianu " << endl;
	//cout << odwrJ[0][0] << "  " << odwrJ[0][1] << endl << odwrJ[1][0] << "   " << odwrJ[1][1] << endl << endl;;

	//WYZNACZNIK JAKOBIANU
	double detJ;
	detJ = (jakobian[0][0] * jakobian[1][1]) - (jakobian[0][1] * jakobian[1][0]);

	//cout << endl << "detJ= " << detJ << endl;
	
	return detJ;
}

void createLocalC(double detJ, double** C, int pkt) {
	double NNT[4][4];
	//N*N(T)  
	for (int i = 0; i < 4; i++)
	{
		NNT[i][0] = el4->N[pkt][i] * el4->N[pkt][0];
		NNT[i][1] = el4->N[pkt][i] * el4->N[pkt][1];
		NNT[i][2] = el4->N[pkt][i] * el4->N[pkt][2];
		NNT[i][3] = el4->N[pkt][i] * el4->N[pkt][3];
		//cout << NNT[i][0] << "    " << NNT[i][1] << "   " << NNT[i][2] << "    " << NNT[i][3] << endl;
	}//cout << endl;
	//cout << "C lokalna: " << endl;
	
	for (int i = 0; i < 4; i++)
	{
		C[i][0] = Dt->c * NNT[i][0] * Dt->ro * detJ;
		C[i][1] = Dt->c * NNT[i][1] * Dt->ro * detJ;
		C[i][2] = Dt->c * NNT[i][2] * Dt->ro * detJ;
		C[i][3] = Dt->c * NNT[i][3] * Dt->ro * detJ;
		//cout << C[i][0] << " " << C[i][1] << " " << C[i][2] << " " << C[i][3] << endl;
			
	}
}

void createHbc(double** Hbc, int numerElementu, int nrPowierzchni){
	//[Hbc]= | alfa*({N}*{N}T)*detHbc
	double NNT[4][4][schematCalkowania];
	for (int pkt = 0; pkt < schematCalkowania; pkt++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++) {
				//dla 0 pow: x= 0,1,2  ||  1 pow: 3,4,5
				//dla 2 pow: x= 6,7,8  ||  3 pow: 9,10,11
				int x = pkt + schematCalkowania * nrPowierzchni; //leci po punktach ca³kowania dla 1 el.skoñczonego
				NNT[i][j][pkt] = el4->Nbc[x][i] * el4->Nbc[x][j] * el4->wagaHbc[pkt];
				Hbc[i][j] = Hbc[i][j]+NNT[i][j][pkt];
				//	cout << NNT[i][j][pc] << "   ";        
			}
		}
	}

	int x = 0;
	if (nrPowierzchni == 0) x = 1;
	if (nrPowierzchni == 1) x = 2;
	if (nrPowierzchni == 2) x = 3;
	if (nrPowierzchni == 3) x = 0;
	double detHbc = sqrt(pow((Nd[Elem[numerElementu].ID[x]].x - Nd[Elem[numerElementu].ID[nrPowierzchni]].x), 2) + 
		pow((Nd[Elem[numerElementu].ID[x]].y - Nd[Elem[numerElementu].ID[nrPowierzchni]].y), 2)) / 2.0;

	for (int i = 0; i < 4; i++)
	{
			Hbc[i][0] = Hbc[i][0]* Dt->alfa * detHbc;
			Hbc[i][1] = Hbc[i][1]* Dt->alfa * detHbc;
			Hbc[i][2] = Hbc[i][2]* Dt->alfa * detHbc;
			Hbc[i][3] = Hbc[i][3]* Dt->alfa * detHbc;
	}
}

void createLocalH(double detJ, double** H, int pkt, double** odwrJ)
{
	double dNdX[4], dNdY[4];
	double hPoX[4][4], hPoY[4][4];

	dNdX[0] = (el4->dNdKsi[pkt][0] * odwrJ[0][0] + el4->dNdEta[pkt][0] * odwrJ[0][1]) * (1 / detJ);
	dNdX[1] = (el4->dNdKsi[pkt][1] * odwrJ[0][0] + el4->dNdEta[pkt][1] * odwrJ[0][1]) * (1 / detJ);
	dNdX[2] = (el4->dNdKsi[pkt][2] * odwrJ[0][0] + el4->dNdEta[pkt][2] * odwrJ[0][1]) * (1 / detJ);
	dNdX[3] = (el4->dNdKsi[pkt][3] * odwrJ[0][0] + el4->dNdEta[pkt][3] * odwrJ[0][1]) * (1 / detJ);

	dNdY[0] = (el4->dNdKsi[pkt][0] * odwrJ[1][0] + el4->dNdEta[pkt][0] * odwrJ[1][1]) * (1 / detJ);
	dNdY[1] = (el4->dNdKsi[pkt][1] * odwrJ[1][0] + el4->dNdEta[pkt][1] * odwrJ[1][1]) * (1 / detJ);
	dNdY[2] = (el4->dNdKsi[pkt][2] * odwrJ[1][0] + el4->dNdEta[pkt][2] * odwrJ[1][1]) * (1 / detJ);
	dNdY[3] = (el4->dNdKsi[pkt][3] * odwrJ[1][0] + el4->dNdEta[pkt][3] * odwrJ[1][1]) * (1 / detJ);

	//dNdX*dNdX(T)  dNdY*dNdY(T)
	for (int i = 0; i < 4; i++)
	{
		hPoX[i][0] = dNdX[i] * dNdX[0];
		hPoX[i][1] = dNdX[i] * dNdX[1];
		hPoX[i][2] = dNdX[i] * dNdX[2];
		hPoX[i][3] = dNdX[i] * dNdX[3];

		hPoY[i][0] = dNdY[i] * dNdY[0];
		hPoY[i][1] = dNdY[i] * dNdY[1];
		hPoY[i][2] = dNdY[i] * dNdY[2];
		hPoY[i][3] = dNdY[i] * dNdY[3];
	}

	for (int i = 0; i < 4; i++)
	{
		H[i][0] = (hPoY[i][0] + hPoX[i][0]) * Dt->k * detJ;
		H[i][1] = (hPoY[i][1] + hPoX[i][1]) * Dt->k * detJ;
		H[i][2] = (hPoY[i][2] + hPoX[i][2]) * Dt->k * detJ;
		H[i][3] = (hPoY[i][3] + hPoX[i][3]) * Dt->k * detJ;
	}
}

void createP(double* P, int numerElementu, int nrPowierzchni) {
	// {P} = | alfa * {N} * temp_otoczenia * detP
	for (int pkt = 0; pkt < schematCalkowania; pkt++)
	{
		for (int i = 0; i < 4; i++)
		{
			//dla 0 pow: x= 0,1,2  ||  1 pow: 3,4,5
			//dla 2 pow: x= 6,7,8  ||  3 pow: 9,10,11
			int x = pkt + schematCalkowania * nrPowierzchni; //leci po punktach ca³kowania
			P[i] += el4->Nbc[x][i] * el4->wagaHbc[pkt];
		}
	}

	int x = 0;
	if (nrPowierzchni == 0) x = 1;
	if (nrPowierzchni == 1) x = 2;
	if (nrPowierzchni == 2) x = 3;
	if (nrPowierzchni == 3) x = 0;

	double detP = sqrt(pow((Nd[Elem[numerElementu].ID[x]].x - Nd[Elem[numerElementu].ID[nrPowierzchni]].x), 2) + pow((Nd[Elem[numerElementu].ID[x]].y - Nd[Elem[numerElementu].ID[nrPowierzchni]].y), 2)) / 2.0;
	
	P[0] = P[0]* Dt->alfa * detP * Dt->Temp_otoczenia;
	P[1] = P[1]* Dt->alfa * detP * Dt->Temp_otoczenia;
	P[2] = P[2]* Dt->alfa * detP * Dt->Temp_otoczenia;
	P[3] = P[3]* Dt->alfa * detP * Dt->Temp_otoczenia;
}

void gauss(double* temp, double** Hglobal, double* Pglobal)
{
	int x, y, z, w_1, w_2;
	double max = 0;
	for (x = 0; x < Dt->nN; x++){
		while (Hglobal[x][x] == 0){
			w_1 = x;
			for (y = 0; y < Dt->nN; y++){
				while (abs(Hglobal[y][x]) > max){
					max = Hglobal[y][x];
					w_2 = y;
				}
			}
		}
		while (max != 0){
			for (z = 0; z < Dt->nN; z++)	{
				swap(Hglobal[w_1][z], Hglobal[w_2][z]);
			} swap(Pglobal[w_1], Pglobal[w_2]);
			max = 0;
		}
	}
	for (z = 0; z < Dt->nN - 1; z++){
		for (x = z + 1; x < Dt->nN; x++){
			Hglobal[x][z] = Hglobal[x][z]/Hglobal[z][z];
			for (y = z + 1; y < Dt->nN + 1; y++){
				if (y == Dt->nN)
					Pglobal[x] = Pglobal[x]- Hglobal[x][z] * Pglobal[z];
				else Hglobal[x][y] = Hglobal[x][y]- Hglobal[x][z] * Hglobal[z][y];
			} Hglobal[x][z] = 0;
		}
	}

	for (x = Dt->nN - 1; x >= 0; x--){
		temp[x] = Pglobal[x];
		for (y = x + 1; y < Dt->nN; y++){
			if (y == Dt->nN) {
				temp[x] = temp[x] - (Pglobal[x] * temp[y]);
			}
			else {
				temp[x] = temp[x] - (Hglobal[x][y] * temp[y]);
			}
		}temp[x] = temp[x]/ Hglobal[x][x];
	}
}

void max_min(int iteracje, double* tem)
{
	double max = 0, min = 0;

	min = tem[0];
	for (int i = 0; i < Dt->nN; i++)
	{
		if (tem[i] < min) min = tem[i];
	}

	tem[0] = max;
	for (int i = 0; i < Dt->nN; i++)
	{
		if (tem[i] > max) max = tem[i];
	}

	cout << iteracje << "              " << min << "             " << max << endl;
}

void createHglobalCglobal()
{
	//tworzenie macierzy H globalnej 
	double** Hglobal = new double* [Dt->nN];
	for (int i = 0; i < Dt->nN; i++)
	{
		Hglobal[i] = new double[Dt->nN];
	}
	//tworzenie macierzy C globalnej 
	double** Cglobal = new double* [Dt->nN];
	for (int i = 0; i < Dt->nN; i++)
	{
		Cglobal[i] = new double[Dt->nN];
	}
	//tworzenie wektroa P globalnego
	double* Pglobal = new double[Dt->nN];

	el4 = new Derivatives();
	double x[4], y[4];
	double w1[size], w2[size];
	double** jakobian = new double* [2];
	jakobian[0] = new double[2];
	jakobian[1] = new double[2];
	double** odwrJ = new double* [2];
	odwrJ[0] = new double[2];
	odwrJ[1] = new double[2];


	double** Hlpc = new double* [4];
	double** Clpc = new double* [4];
	double** Hbc = new double* [4];
	double* P = new double[4];
	Hlpc[0] = new double[4]; Clpc[0] = new double[4]; Hbc[0] = new double[4];
	Hlpc[1] = new double[4]; Clpc[1] = new double[4]; Hbc[1] = new double[4];
	Hlpc[2] = new double[4]; Clpc[2] = new double[4]; Hbc[2] = new double[4];
	Hlpc[3] = new double[4]; Clpc[3] = new double[4]; Hbc[3] = new double[4];

	double Hlocal[4][4];
	double Clocal[4][4];
	double Plocal[4];

	cout << "Time[s]        MinTemp[s]      MaxTemp[s]" << endl;
	for (int iteracje = 0; iteracje < (Dt->timeS / Dt->deltaTau); iteracje++)
  {
	//cout << "Iteracja : " << (iteracje * Dt->deltaTau)+Dt->deltaTau << endl;
	int It = (iteracje * Dt->deltaTau) + Dt->deltaTau;

	//ZEROWANHIE
	for (int i = 0; i < Dt->nN; i++)
	{
		for (int j = 0; j < Dt->nN; j++)
		{
			Hglobal[i][j] = 0;
		}
	}
	for (int i = 0; i < Dt->nN; i++)
	{
		for (int j = 0; j < Dt->nN; j++)
		{
			Cglobal[i][j] = 0;
		}
	}
	for (int i = 0; i < Dt->nN; i++)
	{
		Pglobal[i] = 0;
	}

	//H GLOBALNA + Hbc
	for (int i = 0; i < Dt->nE; i++)
	{
		//zerowanie H lokalnej 
		for (int k = 0; k < 4; k++) {
			for (int j = 0; j < 4; j++) {
				Hlocal[k][j] = 0;
				Hbc[k][j] = 0;
			}
		}

		//POBIERAMY WSPÓ£RZÊDNE WÊZ£ÓW
		x[0] = Nd[Elem[i].ID[0]].x;    y[0] = Nd[Elem[i].ID[0]].y;
		x[1] = Nd[Elem[i].ID[1]].x;    y[1] = Nd[Elem[i].ID[1]].y;
		x[2] = Nd[Elem[i].ID[2]].x;    y[2] = Nd[Elem[i].ID[2]].y;
		x[3] = Nd[Elem[i].ID[3]].x;    y[3] = Nd[Elem[i].ID[3]].y;

		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				for (int z = 0; z < size; z++) { //z- leci po punktach ca³kowania
					double detJ = createJacobian(jakobian, x, y, z, odwrJ);
					createLocalH(detJ, Hlpc, z, odwrJ);
					Hlocal[j][k] += Hlpc[j][k] * el4->w1[z] * el4->w2[z];

				}
			}
		}
		for (int nrPow = 0; nrPow < 4; nrPow++) {  
			//cout << nrPow << endl;
			createHbc(Hbc, i, nrPow);  //Hbc, nr Elem, nrPowierzchni,
			int x = 0;
			if (nrPow== 0) x = 1;
			if (nrPow == 1) x = 2;
			if (nrPow == 2) x = 3;
			if (nrPow == 3) x = 0;
			for (int a = 0; a < 4; a++)
			{
				for (int b = 0; b < 4; b++)
				{
					if (Nd[Elem[i].ID[nrPow]].Bc == 1 && Nd[Elem[i].ID[x]].Bc == 1) {
						Hlocal[a][b] += Hbc[a][b];
						//cout << "Warunke spe³niony"<<endl;
					}

					   Hbc[a][b] = 0.0;
				}

			}
		}

		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				Elem[i].H[j][k] = Hlocal[j][k];
			}
		}
		/*cout << "Numer elementu: " << i << endl;
		for (int i = 0; i < 4; i++) {
			cout << Hlocal[0][i] << " " << Hlocal[1][i] << "  " << Hlocal[2][i] << "  " << Hlocal[3][i] << endl;
		}cout << endl;*/
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				//dodawanie Hlokalnych w miejsca wspolrzednych globalnych wezlow tych elem.
				Hglobal[Elem[i].ID[j]][Elem[i].ID[k]] += Elem[i].H[j][k];

			}
		}
	}

	//C GLOBALNA 
	for (int i = 0; i < Dt->nE; i++)
	{
		//zerowanie C lokalnej 
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Clocal[i][j] = 0;
			}
		}

		//POBIERAMY WSPÓ£RZÊDNE WÊZ£ÓW
		x[0] = Nd[Elem[i].ID[0]].x;    y[0] = Nd[Elem[i].ID[0]].y;
		x[1] = Nd[Elem[i].ID[1]].x;    y[1] = Nd[Elem[i].ID[1]].y;
		x[2] = Nd[Elem[i].ID[2]].x;    y[2] = Nd[Elem[i].ID[2]].y;
		x[3] = Nd[Elem[i].ID[3]].x;    y[3] = Nd[Elem[i].ID[3]].y;

		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				for (int z = 0; z < size; z++) {
					double detJ = createJacobian(jakobian, x, y, z, odwrJ);
					createLocalC(detJ, Clpc, z);
					Clocal[j][k] += Clpc[j][k]*el4->w1[z]*el4->w2[z];
					Elem[i].C[j][k] = Clocal[j][k];
				}
			}
		}

		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				Cglobal[Elem[i].ID[j]][Elem[i].ID[k]] += Elem[i].C[j][k];
			}
		}
	}

	//oblicz P global
	for (int i = 0; i < Dt->nE; i++)
	{
		//zerowanie P lokalnej 
		for (int j = 0; j < 4; j++) {
			Plocal[j] = 0;
			P[j] = 0;
		}

		for (int nrPow = 0; nrPow < 4; nrPow++) {
			createP(P, i, nrPow);
			for (int a = 0; a < 4; a++)
			{
				int x = 0;
				if (nrPow == 0) x = 1;
				if (nrPow == 1) x = 2;
				if (nrPow == 2) x = 3;
				if (nrPow == 3) x = 0;
				if (Nd[Elem[i].ID[nrPow]].Bc == 1 && Nd[Elem[i].ID[x]].Bc == 1)
					Plocal[a] += P[a];
				P[a] = 0;
			}
		}

		for (int j = 0; j < 4; j++) {
			Elem[i].P[j] = Plocal[j];
		}

		for (int j = 0; j < 4; j++)
		{
			//wybieram z elementu konkretne ID i dopisuje to co siedzi w elemencie Plocalnym na miejscu
			Pglobal[Elem[i].ID[j]] += Elem[i].P[j];
		}
	}


	//wypisywanie H, C, P globalnych
	/*
	cout << "HGlobalna: " << endl;
	for (int i = 0; i < Dt->nN; i++)
	{
		for (int j = 0; j < Dt->nN; j++)
		{
			cout << Hglobal[i][j] << " ";
		} cout << endl;
	}
	
	cout << "PGlobalna: " << endl;
	for (int i = 0; i < Dt->nN; i++)
	{
		cout << Pglobal[i] << " ";

	}

	cout << endl << endl;
	cout << "Clobalna: " << endl;
	for (int i = 0; i < Dt->nN; i++)
	{
		for (int j = 0; j < Dt->nN; j++)
		{
			cout << Cglobal[i][j] << " ";
		} cout << endl;
	}
	
	*/


	
	//ZASTÊPCZE  H
	// [C]/dT
	for (int i = 0; i < Dt->nN; i++) // -->ilosci wêz³ów
	{
		for (int j = 0; j < Dt->nN; j++)
		{
			Cglobal[i][j] = Cglobal[i][j] / Dt->deltaTau;
		}
	}
	
	//[H]=[H]+[C]/dT
	for (int i = 0; i < Dt->nN; i++)
	{
		for (int j = 0; j < Dt->nN; j++)
		{
			Hglobal[i][j] = Hglobal[i][j] + Cglobal[i][j];
			//cout << Hglobal[i][j] << "  ";
		}
		//cout << endl;
	}
	/*
	cout << "H zastepcze ([H]=[H]+[C]/dT): " << endl;
	for (int i = 0; i < Dt->nN; i++)
	{
		for (int j = 0; j < Dt->nN; j++)
		{
			cout << Hglobal[i][j] << " ";
		} cout << endl;
	}
	*/
	//ZASTÊPCZE  P
	double* CDTT0 = new double[Dt->nN];

	for (int i = 0; i < Dt->nN; i++)
	{
		CDTT0[i] = 0;
		for (int j = 0; j < Dt->nN; j++) {
			//  ([C] / dT) * t0
			CDTT0[i] = CDTT0[i]+Cglobal[i][j] * Nd[j].t0; //ka¿dy wêze³ mno¿ymy przez t0
		}
	} 

	//{P} = {P} + ([C] / dT) * t0
	for (int i = 0; i < Dt->nN; i++)
	{
		Pglobal[i] = Pglobal[i] + CDTT0[i];
	}

	/*
	cout << "P zastepcze ({P} = {P} + ([C] / dT) * t0) : " << endl;
	for (int i = 0; i < Dt->nN; i++)
	{
		cout << Pglobal[i] << " ";

	}cout << endl;
	*/
	//ROZWI¥ZANIE UK£ADU RÓWNAÑ
	// [H]*{t}=-{P}
	double* temp = new double[Dt->nN];
	gauss(temp,Hglobal,Pglobal); //tu oblicza nowe temp w wezlach 
	for (int i = 0; i < Dt->nN; i++)
	{
		Nd[i].t0 = temp[i]; //tu przypisujemy nowe temp (nowe staja sie starymi w nowym kroku)
	}

     max_min(It, temp);
}

}

int main() {

	creatingMesh();
	createHglobalCglobal();

	system("PAUSE");
	return 0;
}





