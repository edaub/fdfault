#ifndef OUTPUTUNITCLASSHEADERDEF
#define OUTPUTUNITCLASSHEADERDEF

class outputunit
{
public:
	outputunit();
	~outputunit();
    outputunit* get_next_unit() const ;
    void set_next_unit(outputunit* nextunit);
    void write_unit() const;
private:
    double* f;
    outputunit* next;
};

#endif