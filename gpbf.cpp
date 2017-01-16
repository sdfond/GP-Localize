//#include "gpbf.h"
#include "foxDieter.h"

GPBF gpbf;
int main(int argc, char *argv[])
{
	if (argc == 0)
	{
		printf("please specify config name!\n");
		return 0;
	}
	gpbf.parse(argv[1]);
	printf("parsing input.. done.\n");
	gpbf.init();
	printf("initializing.. done.\n");
	gpbf.readData();
	printf("read data.. done.\n");
	gpbf.getSupSet();
	printf("generate support set.. done\n");
	gpbf.filter();
	
	return 0;
}
