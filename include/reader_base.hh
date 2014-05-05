#ifndef __READERCLASS__
#define __READERCLASS__

class cParser{
    public:
	enum readerType {
	    d2type
	};
	virtual void Run() = 0;
	cParser();
	~cParser();
	cParser(readerType type);
        void Print();
    private:
	readerType rtype;
};

#endif
