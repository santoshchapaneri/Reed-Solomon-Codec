typedef unsigned short GFElemType;
typedef unsigned short GFPowRep;

unsigned int FIELD;
GFPowRep POWER;
GFPowRep MSG;
GFPowRep PAR;

struct Field
{
	unsigned   size;
	GFElemType *elements;
	GFPowRep   *powrep;
};

