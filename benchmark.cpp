#include	<sys/time.h>
#include	<iostream>
#include	"wynn.hpp"
using	namespace	wynn;

template<uint64_t	N,	uint64_t	M>
void	benchmark_tiger(void){
	Data<N>	w,g,m;	w.rand();	g.rand();	m.rand();
	timeval	beg,	end;
	gettimeofday(&beg,NULL);
	for(size_t	i=0;	i<M;	i++)	tiger<N>(w.data,g.data,m.data,0.001f);
	gettimeofday(&end,NULL);
	double	t=end.tv_sec-beg.tv_sec+1e-6*(end.tv_usec-beg.tv_usec);
	cerr<<"Tiger\t"<<N/t/1024/1024/1024*M<<" G\n";
}

template<uint64_t	N,	uint64_t	M>
void	benchmark_linear(void){
	Data<N*N>	inp,out,gra;	linear<N,N,N>	w;	inp.rand();
	timeval	beg,	end;
	gettimeofday(&beg,NULL);
	for(size_t	i=0;	i<M;	i++){	w.fw(inp,out);	w.bk(inp,out,gra);	}
	gettimeofday(&end,NULL);
	double	t=end.tv_sec-beg.tv_sec+1e-6*(end.tv_usec-beg.tv_usec);
	cerr<<"linear\t"<<6.0*N*N*N/t/1024/1024/1024*M<<" G\n";
}

int	main(int	ac,	char	**av){
	benchmark_tiger<0x10000,0x100000>();
	benchmark_linear<0x200,0x1000>();
	return	0;
}
