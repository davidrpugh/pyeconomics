% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
atypes = {'lin','pwcub'};
aletter = 'LC';
nmoms = [1 2];
its = [0 1];
ds = [0 1];
js = [1 2 3];


FF = fopen('res/euleraggr.txt','w');
GG = fopen('res/accuracy.txt','w');

textiter = {'Dstst','Dsimul'}

setrefdist;
for i=1:length(nmoms);
  nMoments = nmoms(i);
  for j=2:length(atypes);
    approxtype = atypes{j};
    for l=1:length(its);
      iter=its(l);
      descr = sprintf('Mom%d %s  ',nMoments,textiter{iter+1});
      fprintf(FF,'%s',descr);
      fprintf(GG,'%s',descr);
      for lj=1:length(js);
	whichShocks = js(lj);
	for m=1:length(ds);
	  useDsel = ds(m);
	  tagSim = sprintf('%s%d%d%d_%d',...
			   aletter(j),whichShocks,nMoments,iter,useDsel);
	  cmd = sprintf('load res/ser%s; ser%d=ser%s;id%d=id%s; ',...
			tagSim,useDsel,tagSim,useDsel,tagSim);
	  eval(cmd);
	end
	%indx = 1:min(100000,size(ser0,1));
	indx = 1:size(ser0,1);
	fprintf(FF,' & %14s',mynumb(mean(ser0(indx,id0.maxres))),mynumb(mean(ser0(indx,id0.meanres))));
	mdiff = abs(ser0(indx,1)./ser1(indx,1)-1);
	fprintf(GG,' & %14s',mynumb(max(mdiff)),mynumb(mean(mdiff)));

	tagSimSh = sprintf('%s%d%d%d',...
			   aletter(j),whichShocks,nMoments,iter);
	HH = fopen(sprintf('res/acc%s',tagSimSh),'w');
	for t=1:size(ser1,1)
	  fprintf(HH,'%d  %18.10f  %18.10f  %18.10f  %18.10f\n',...
		  ser0(t,id1.ia),ser0(t,id1.meanU),ser0(t,id1.meanE),ser1(t,id1.meanU),ser1(t,id1.meanE));
	end
	fclose(HH);
	JJ = fopen(sprintf('res/tsacc%s',tagSimSh),'w');
      end
      fprintf(FF,'\\\\\n');
      fprintf(GG,'\\\\\n');
    end
  end
end


fclose(FF);
fclose(GG);
