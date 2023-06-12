#include	<sys/time.h>
#include	<iostream>
#include	"wynn.hpp"
using	namespace	wynn;

template<unsigned	N,	unsigned	M>
void	benchmark_tiger(void){
	float	*w=new	float[N],*g=new	float[N],*m=new	float[N];
	timeval	beg,	end;
	gettimeofday(&beg,NULL);
	for(size_t	i=0;	i<M;	i++)	tiger<N>(w,g,m,0.001f);
	gettimeofday(&end,NULL);
	double	t=end.tv_sec-beg.tv_sec+1e-6*(end.tv_usec-beg.tv_usec);
	cerr<<"Tiger\t"<<N/t/1024/1024/1024*M<<" G\n";
	delete	[]	w;	delete	[]	g;	delete	[]	m;
}

int	main(int	ac,	char	**av){
	benchmark_tiger<0x100000,0x100000>();
	return	0;
}
