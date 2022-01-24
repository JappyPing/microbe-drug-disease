clc;tic;
load interaction.mat;
[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number

circSim01  = GSM( interaction'); % microbe gaussian
disSim01  = GSD( interaction' ); % disease gaussian

disSim02  = cosSim( interaction); % disease cosine
circSim02  = cosSim( interaction');   % microbe cosine

w1 = 0.6;
w2 = 0.3;
circRNA_sim  = LKF(w1,circSim01,circSim02); % LKF microbe gaussian and cosine
dis_sim  = LKF(w2,disSim01,disSim02); % LKF disease gaussian and cosine 

F = NCPLDA(circRNA_sim,dis_sim,interaction');

%%
microbe_list = importdata('microbe_list.txt','%s');
disease_list = importdata('disease_list.txt');

%%
[F,LP_rank_known] =Rank( F, interaction',microbe_list,disease_list);
Write_file( F );
toc;
