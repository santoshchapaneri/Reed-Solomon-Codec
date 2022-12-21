
//#define DEBUG

void ZechsLogarithms();
GFPowRep GFMultiplication(GFPowRep elem1, GFPowRep elem2);
GFPowRep GFAddition(GFPowRep elem1, GFPowRep elem2);
GFPowRep GFExponentiation(GFPowRep elem1, int power);
void GeneratorCreator(GFPowRep *GenPoly);

void RSEncoder(char *fin, char *fout);
void RSDecoder(char *fin, char *fout);

GFPowRep GFPolyEvaluate(GFPowRep *Poly, int degree, GFPowRep arg);

