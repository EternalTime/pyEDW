%   Parameter setup
%   theta (Parameter Vector)

f       =   0.88;
Ag      =   0.3;
Ab      =   0.1;
Aw      =   0.6;
Topt    =   300;
Q       =   0.1;
gG      =   1;
gD      =   0.2;
tauS    =   3;
tauE    =   5;
dT      =   30;

dt = 0.1;

theta0 = [ ...
            f , ...             %    1:     f
            1-Ag , ...          %    2:     1 - A_G
            Ab-Ag, ...          %    3:     A_B - A_G
            Aw-Ag, ...          %    4:     A_W - A_G
            Q, ...              %    5:     Q
            gD/gG, ...          %    6:     gamma_D / gamma_G
            1/(gG*tauE), ...    %    7:     1 / gamma_G tau_E
            1/(gG*tauS), ...    %    8:     1 / gamma_G tau_S
            8*(Topt/dT)^4, ...  %    9:     8 (T_opt / Delta T)^4
            0.05, ...           %    10:    delta
            0 ...               %    11:    lambda
        ];

dTs = linspace(0,80,129);
dTs(1) = [];
lambda = linspace(-0.7,1.4,400);
Ls = (1+lambda);
Ts = Ls.^0.25;
N = 500;

pc = parcluster('local')
JOB_ID = getenv('SLURM_JOBID')
CPUS = getenv('SLURM_CPUS_PER_TASK')
pc.JobStorageLocation = strcat('/local_scratch/',JOB_ID)
parpool(pc,str2num(CPUS))

for ii = str2num(getenv('SLURM_ARRAY_TASK_ID'))

    theta = theta0;
    theta(9) = 8*(Topt/dTs(ii))^4;

    data = zeros([length(Ls),N,4]);
    t1 = tic;
    for ll = 1:length(lambda)
        tic

        theta(11) = lambda(ll);

        for nn = 1:N
            f1 = rand()*f; f2 = 1 - f1;
            y = rand();
            f1 = y*f1; f2 = y*f2;

            x = [f1 f2 Ts(ll) Ls(ll)];
            for tt = 1:2000
                x = updateExoDaisyWorld(x,dt,theta);
            end
            data(ll,nn,:) = x;
        end
        save("Data/data_"+num2str(ii,"%03.0f")+".mat","data","theta","lambda","Ls","Ts")
        disp(num2str(toc-t1,"%.2f"))
    end
    save("Data/data_"+num2str(ii,"%03.0f")+".mat","data","theta","lambda","Ls","Ts")
end

delete(gcp('nocreate'))