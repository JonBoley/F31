#include "model.hpp"
int main(int argc,char* argv[])
{
   //This automatically calls the function SingleAN(argc,argv) before run()
   SingleAN *test = new SingleAN(argc,argv);
   test->run();
};
