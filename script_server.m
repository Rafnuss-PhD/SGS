%% SECTION: COMPUTATIONAL SAVING

clear all;

% General setting
parm.k.covar(1).model = 'spherical';
parm.k.covar(1).azimuth = 0;
parm.k.covar(1).c0 = 1;
parm.k.covar(1).alpha = 1;
parm.seed_path = 'default';
parm.seed_search = 'shuffle';
parm.seed_U = 'default';

parm.k.covar(1).range0 = [15 15];
parm.k.wradius = 3;
parm.n_real = 1;
parm.saveit = 0;

parm.mg = 1;


N=2.^[8 10 13]+1;%[255 511 1023 2047 4095 8191];
K=[108];

%nt=[50 50 50];

for i_n=1:numel(N)
    for i_k=1:numel(K)
        nx = N(i_n); % no multigrid
        ny = N(i_n);
        parm.k.nb = K(i_k);
        
        %tt={};tt.global=nan(nt(i_n),1);tt.real=nan(nt(i_n),1);
        %for i=1:nt(i_n)
            [~,t] = SGS_cst_par(nx,ny,parm);
         %   tt.global(i) = t.global;
         %   tt.real(i) = t.real;
        %end
        title = ['T_2_cst_par_' num2str(N(i_n)) 'N_' num2str(K(i_k)) 'K' ];
        save(['./cst_path_paper/' title ],'parm','t')
        
        %tt={};tt.global=nan(nt(i_n),1);
        %for i=1:nt(i_n)
            [~,t] = SGS_trad(nx,ny,parm);
        %    tt.global(i) = t.global;
        %   tt.real(i) = t.real;
        %end
        
        title = ['T_2_trad_' num2str(N(i_n)) 'N_' num2str(K(i_k)) 'K' ];
        save(['./cst_path_paper/' title ],'parm','t')
    end
end






