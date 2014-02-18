function [bigevent] = bigevents(amplitudethreshold,timelimit,tracespikedel,time,varargin)
% Calculate the big events in a trace.
%   INPUT: 
%       amplitudethreshold: compute events only larger than that.
%                           (typically 2-3 mV in order not to blow up the
%                           calculations!)
%       timelimit:          how fast the event has to be. typically 4 ms
%                           (smaller than the membrane tau)
%       tracespikedel:      the trace with the spikes already chopped
%       time:               just a vector containing the time of each step
%                           (same size of tracespikedel)
%       
%   OUTPUT: 
%           bigevent.ampli contains the amplitude of the events
%           bigevent.binstart contains the bin of the onset of big events
%           bigevent.binduration contains the duration of the big events
%           bigevent.voltstart contains the starting voltage of the event
%           

TIMEPERBIN=time(2)-time(1);

binlimit=round(timelimit./TIMEPERBIN);

deltav=tracespikedel(binlimit:end)-tracespikedel(1:end-binlimit+1);

#amplitudethreshold
#timelimit
#tracespikedel(1)
#time(2)

suprathrev=find(deltav>amplitudethreshold);

diffsuprathrev=diff(suprathrev);
blockstart=find(diffsuprathrev>1);
blockstart=[1; blockstart+1];
%blockstart(end)=[];

timemin=0.005;
binmin=round(timemin./TIMEPERBIN);

difftime=time(suprathrev(blockstart(2:end)))-time(suprathrev(blockstart(1:end-1)));
closestarts=find(difftime<timemin);
blockstart(closestarts+1)=[];

blockend=blockstart(2:end)-1;
blockend=[blockend; length(suprathrev)];

nevents=length(blockstart);
ii=0;

while ii<nevents
   ii=ii+1;
   
   npoints=blockend(ii)-blockstart(ii)+1;
   jj=0;
   maxvoltdiff=[];stopindex=[];
   while jj<npoints
       jj=jj+1;
       
       stb_a=suprathrev(blockstart(ii))+jj;
       stb_b=min(suprathrev(blockstart(ii))+jj+binlimit,length(tracespikedel));
       [maxvoltdiff(jj),stopindex(jj)]=max(tracespikedel( stb_a  :  stb_b)-tracespikedel( stb_a  ));
   end
   
   if jj==0
       jj=1;
   end
   
   [bigevent.ampli(ii),startindex]=max(maxvoltdiff);
   bigevent.binstart(ii)=suprathrev(blockstart(ii))+startindex;
   bigevent.binduration(ii)=stopindex(startindex)-1;
   bigevent.voltstart(ii)=tracespikedel(suprathrev(blockstart(ii))+startindex);
   
   binstop=bigevent.binstart(ii)+bigevent.binduration(ii);  
   
end

end

