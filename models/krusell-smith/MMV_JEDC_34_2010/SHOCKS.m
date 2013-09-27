% The program for the article "Solving the incomplete markets model with
% aggregate uncertainty using the Krusell-Smith algorithm" from the special 
% JEDC issue edited by Den Haan, Judd and Juillard (2008)  
%
% Written by Lilia Maliar, Serguei Maliar and Fernando Valli (2008)

function [idshock,agshock]  = SHOCKS(prob,T,N,ur_b);

disp('Generating shocks');

idshock=zeros(T,N); % matrix of idiosyncratic shocks 
agshock=zeros(T,1); % vector of aggregate shocks

%__________________________________________________________________________
%
% Transition probabilities between the aggregate states 
%__________________________________________________________________________

% prob_ag(i,j) is the probability of tomorrow's agg. shock (i=1,2) given 
% today's agg. shock (j=1,2)

prob_ag=zeros(2,2);  
prob_ag(1,1)=prob(1,1)+prob(1,2); prob_ag(2,1)=1-prob_ag(1,1);  
prob_ag(2,2)=prob(3,3)+prob(3,4); prob_ag(1,2)=1-prob_ag(2,2);

%__________________________________________________________________________
%
% Probability of an idiosyncratic shock epsilon' given that aggregate shock
% s' is realized;
%__________________________________________________________________________

% p_xy_zw is the probability of idiosyncratic shock epsilon'=w conditional 
% on aggregate shocks s'=y, s=x and idiosyncratic shock epsilon=z

p_bb_uu = prob(1,1)/prob_ag(1,1); p_bb_ue=1-p_bb_uu;
p_bb_ee = prob(2,2)/prob_ag(1,1); p_bb_eu=1-p_bb_ee;
p_bg_uu = prob(1,3)/prob_ag(2,1); p_bg_ue=1-p_bg_uu;
p_bg_ee = prob(2,4)/prob_ag(2,1); p_bg_eu=1-p_bg_ee;
p_gb_uu = prob(3,1)/prob_ag(1,2); p_gb_ue=1-p_gb_uu;
p_gb_ee = prob(4,2)/prob_ag(1,2); p_gb_eu=1-p_gb_ee;
p_gg_uu = prob(3,3)/prob_ag(2,2); p_gg_ue=1-p_gg_uu;
p_gg_ee = prob(4,4)/prob_ag(2,2); p_gg_eu=1-p_gg_ee;
%__________________________________________________________________________
%
% Generation of the aggregate shocks 
%__________________________________________________________________________

agshock(1)=1; % assume that initially (t=1), the economy is in a bad state 

% To generate shocks in subsequent periods, we draw random numbers. If a 
% random number drawn is <= probability of a bad shock conditional on 
% agshock(t-1)), then set agshock(t)=1; otherwise set agshock(t)=2
for t=2:T
   raux=rand; 
   if raux<=prob_ag(1,agshock(t-1)) 
      agshock(t)=1; % 
   else
      agshock(t)=2;
   end
end

%__________________________________________________________________________
%
% Generation of the idiosyncratic shocks for all agents in the first period
%__________________________________________________________________________

for i=1:N
   raux=rand;
   if raux<=ur_b % in the first period agg. shock is bad; if a  random number 
       % drawn is <= the probability of being unemployed in a bad agg. state, 
       % then set idshock(t)=1; otherwise set idshock(t)=2
      idshock(1,i)=1;
   else
      idshock(1,i)=2;
   end
end
%__________________________________________________________________________
%
% Generation of the idiosyncratic shocks for all agents starting from the 
% second period
%__________________________________________________________________________

for t=2:T
      
   if agshock(t-1)==1 & agshock(t)==1 % if the previous agg. shock was bad 
                                      % and the current agg. shock is bad
      for i=1:N
         raux=rand;
         if idshock(t-1,i)==1 % if the previous idiosyncratic shock was 1 
            if raux<=p_bb_uu
               idshock(t,i)=1;
            else
               idshock(t,i)=2;
            end
         else                 % if the previous idiosyncratic shock was 2 
            if raux<=p_bb_ee
               idshock(t,i)=2;
            else
               idshock(t,i)=1;
            end
         end
      end
   end
   
   if agshock(t-1)==1 & agshock(t)==2 % if the previous agg. shock was bad 
                                      % and the current agg. shock is good
      for i=1:N
         raux=rand;
         if idshock(t-1,i)==1
            if raux<=p_bg_uu
               idshock(t,i)=1;
            else
               idshock(t,i)=2;
            end
         else
            if raux<=p_bg_ee
               idshock(t,i)=2;
            else
               idshock(t,i)=1;
            end
         end
      end
   end
   
   if agshock(t-1)==2 & agshock(t)==1 % if the previous agg. shock was good 
                                      % and the current agg. shock is bad
      for i=1:N
         raux=rand;
         if idshock(t-1,i)==1
            if raux<=p_gb_uu
               idshock(t,i)=1;
            else
               idshock(t,i)=2;
            end
         else
            if raux<=p_gb_ee
               idshock(t,i)=2;
            else
               idshock(t,i)=1;
            end
         end
      end
   end
   
   if agshock(t-1)==2 & agshock(t)==2 % if the previous agg. shock was good 
                                      % and the current agg. shock is good
      for i=1:N
         raux=rand;
         if idshock(t-1,i)==1
            if raux<=p_gg_uu
               idshock(t,i)=1;
            else
               idshock(t,i)=2;
            end
         else
            if raux<=p_gg_ee
               idshock(t,i)=2;
            else
               idshock(t,i)=1;
            end
         end
      end
   end
end