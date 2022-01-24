%% run 10-fold cross validation for methods of BNNR, DR2DI, MBiRW, MSBMF

function [] = ten_fold_CV(num)
    addpath('../Data/')
    %% 1 load data
    addpath('../Data/BNNR')
    load Fdataset
    Wrr = drug;
    Wdd = disease;
    Wdr = didr;
    Wrd = Wdr';
    [nd,nr] = size(Wdr);

    addpath('../Data/DR2DI')
    load Disease_ID
    DiseaseIDs=disease_id;
    load Drug_id
    DrugIDs=drug_id;
    ds=load('drug_structure_Sim.mat');
    ds=cell2mat(struct2cell(ds));
    td=load('target_domain_Sim.mat');
    td=cell2mat(struct2cell(td));
    tg=load('target_go_Sim.mat');
    tg=cell2mat(struct2cell(tg));
    [row,col]=size(ds);
    drug=zeros(row,row);
    for i=1:row
        for j=1:col
            drug(i,j)=(td(i,j)+ds(i,j)+tg(i,j))/3;
        end
    end
    disease=load('disease_Phenotype.mat');
    disease=cell2mat(struct2cell(disease));

    addpath('../Data/MBiRW')
    Wrr1 = importdata('DrugSimMat');
    Wdd1 = importdata('DiseaseSimMat');
    Wrname = importdata('DrugsName');
    Wdname = importdata('DiseasesName');
    newWrr = importdata('shareWrr.mat');
    newWdd = importdata('shareWdd.mat');
    LWrr = Wrr1;
    LWdd = Wdd1;

    addpath('../Data/MSBMF')
    load Fdataset_ms
    Wrr11 = drug_ChemS;
    Wrr22 = drug_AtcS;
    Wrr33 = drug_SideS;
    Wrr44= drug_DDIS;
    Wrr55 = drug_TargetS;
    Wrr2 = [Wrr11, Wrr22, Wrr33, Wrr44, Wrr55];
    Wdd11 = disease_PhS;
    Wdd22 = disease_DoS;
    Wdd2 = [Wdd11, Wdd22];
    min_mn = min(nd, nr);

    known=textread('disease_drug.txt');
    [pp1,qq1]=size(known);
    load unknown_dis_dr
    unknown=a;
    [pp2,qq2]=size(unknown);
    x1=randperm(pp1)';
    x2=randperm(pp2)';

    F1=zeros(nd,nr);
    F2=zeros(nd,nr);
    F3=zeros(nd,nr);
    F4=zeros(nd,nr);
    for cv=1:10
        interaction = Wdr;
        if cv<10
            B1=known(x1((cv-1)*floor(pp1/10)+1:floor(pp1/10)*cv),:);
            B2=unknown(x2((cv-1)*floor(pp2/10)+1:floor(pp2/10)*cv),:);
            for i=1:floor(pp1/10)
                interaction(B1(i,1),B1(i,2))=0;
            end
        else
            B1=known(x1((cv-1)*floor(pp1/10)+1:pp1),:);
            B2=unknown(x2((cv-1)*floor(pp2/10)+1:pp2),:);
            for i=1:pp1-floor(pp1/10)*9
                interaction(B1(i,1),B1(i,2))=0;
            end
        end

        %% BNNR
        maxiter = 300;
        alpha = 1;
        beta = 10;
        tol1 = 2*1e-3;
        tol2 = 1*1e-5;
        T = [Wrr, interaction'; interaction, Wdd];
        [t1, t2] = size(T);
        trIndex = double(T ~= 0);
        [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
        M_recovery = WW((t1-nd+1) : t1, 1 : nr);
        [b1,bb1]=size(B1);
        [b2,bb2]=size(B2);
        for i=1:b1
            F1(B1(i,1),B1(i,2))=M_recovery(B1(i,1),B1(i,2));
        end
        for i=1:b2
            F1(B2(i,1),B2(i,2))=M_recovery(B2(i,1),B2(i,2)); 
        end

        %% DR2DI
        Rt=rls_kron_for_predict(interaction',drug,disease,DiseaseIDs,DrugIDs,1);
        RT=Rt';
        [b1,bb1]=size(B1);
        [b2,bb2]=size(B2);
        for i=1:b1
            F2(B1(i,1),B1(i,2))=RT(B1(i,1),B1(i,2));
        end
        for i=1:b2
            F2(B2(i,1),B2(i,2))=RT(B2(i,1),B2(i,2)); 
        end

        %% MBiRW
        alpha = 0.3;
        l = 2;
        r = 2;
        d = log(9999);
        cr = setparFun(interaction',LWrr);
        cd = setparFun(interaction,LWdd);

        LWrr = 1./(1+exp(cr*LWrr+d));
        LWdd = 1./(1+exp(cd*LWdd+d));

        [RWrr,RWdd] = nManiCluester(LWrr,LWdd,newWrr,newWdd,Wrname,Wdname);

        normWrr = normFun(RWrr);
        normWdd = normFun(RWdd);
        A2=interaction';
        R0 = A2/sum(A2(:));
        Rt1 = R0;

        for t=1:max(l,r)
            ftl = 0;
            ftr = 0;
            
            if(t<=l)
                nRtleft = alpha * normWrr*Rt1 + (1-alpha)*R0;
                ftl = 1;
            end
            if(t<=r)
                nRtright = alpha * Rt1 * normWdd + (1-alpha)*R0;
                ftr = 1;
            end
            Rt1 =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
        end   
        RT1=Rt1';

        [b1,bb1]=size(B1);
        [b2,bb2]=size(B2);
        for i=1:b1
            F3(B1(i,1),B1(i,2))=RT1(B1(i,1),B1(i,2));
        end
        for i=1:b2
            F3(B2(i,1),B2(i,2))=RT1(B2(i,1),B2(i,2)); 
        end

        %% MSBMF
        lambda1 = 0.1;
        lambda2 = 0.01;
        lambda3 = lambda2;
        k = floor(min_mn * 0.7);
        maxiter = 300;
        tol1 = 2*1e-3;
        tol2 = 1*1e-4;
        [U, V, iter] = MSBMF(interaction, Wdd2, Wrr2, lambda1, lambda2, lambda3, k, tol1, tol2, maxiter);
        M_recovery1 = U * V';
        [b1,bb1]=size(B1);
        [b2,bb2]=size(B2);
        for i=1:b1
            F4(B1(i,1),B1(i,2))=M_recovery1(B1(i,1),B1(i,2));
        end
        for i=1:b2
            F4(B2(i,1),B2(i,2))=M_recovery1(B2(i,1),B2(i,2)); 
        end
    end
    str1=strcat('../ten_fold_predict_result/BNNR/Predict_result',num2str(num));
    str2=strcat('../ten_fold_predict_result/DR2DI/Predict_result',num2str(num));
    str3=strcat('../ten_fold_predict_result/MBiRW/Predict_result',num2str(num));
    str4=strcat('../ten_fold_predict_result/MSBMF/Predict_result',num2str(num));
    save(str1,'F1')
    save(str2,'F2')
    save(str3,'F3')
    save(str4,'F4')
end

