function [hw,h_sigma1,h_sigma2,w_max,sigma_x1,sigma_x2,sigma_y1,sigma_y2] = EET_rettangolare_appoggiata_carico_distribuito(a,b,h_1,h_2,t,psi,G,q,alpha,b_x,b_y)

%Restituisce i valori in [mm] degli spessori effettivi per deflessione e
%tensione massime per il primo e secondo strato di vetro
% (risp. "hw" , "h_sigma1" e "h_sigma2"), la deflessione massima al centro della piastra in [mm] ("w_max")
% e le tensione massime in direzione x e y per il primo e secondo strato di vetro in [Mpa] 
% (risp. "sigma_x1" , "sigma_x2" , "sigma_y1" , "sigma_y2").

%Le grandezze di input richieste sono i lati della piastra rettangolare in [mm] ("a" e "b"),
% gli spessori dei due strati di vetro in [mm] ("h_1" e "h_2"), lo spessore dell'intercalare in [mm] ("t"), il 
% coefficiente del grado di accoppiamento tra le lastre di vetro ("psi")
% tabulato, il modulo a taglio dell'intercalare in [Mpa] ("G") ed il carico
% in [N/mm^2] ("q"); sono infine necessari i coefficienti alpha, beta_x e
% beta_y per il calcolo della deflessione e tensioni massime (risp.
% "alpha","b_x" e "b_y")

lambda=b/a;

nu_glass=0.22;
E_glass=70*10^3; %MPa

H=t+((h_1+h_2)/2);
hs_1=H*h_1/(h_1+h_2);
hs_2=H*h_2/(h_1+h_2);
A1=b*h_1;
A2=b*h_2;
Is=(1/b)*(A1*A2*H^2)/(A1+A2);

D_1=E_glass*h_1^3/(12*(1-nu_glass^2));
D_2=E_glass*h_2^3/(12*(1-nu_glass^2));

D_tot=D_1+D_2+(12*D_1*D_2*H^2/((D_1*h_2^2)+(D_2*h_1^2)));

eta=1./(1.+((t./G).*((D_1+D_2)/D_tot)*(12*D_1*D_2*psi)/((D_1*h_2^2)+(D_2*h_1^2))));

D_R=((eta./D_tot)+((1-eta)./(D_1+D_2))).^(-1); %rigidezza effettiva equivalente

hw=(12.*D_R.*(1-nu_glass^2)./E_glass).^(1/3);%spessore effettivo equivalente della lastra monolitica equivalente
h_sigma1=1./((2.*eta*hs_2/(h_1^3+h_2^3+12*Is))+(h_1./hw.^3)).^0.5;
h_sigma2=1./((2.*eta*hs_1/(h_1^3+h_2^3+12*Is))+(h_2./hw.^3)).^0.5;

w_max=alpha*q*(a^4)./D_R; %mm deflessione massima soluzione di Lévy

mx=b_x*q*a^2;%momenti per unità di lunghezza, beta(lambda=1)
my=b_y*q*a^2;

sigma_x1=12*mx*h_sigma1/2./h_sigma1.^3;
sigma_x2=12*mx*h_sigma2/2./h_sigma2.^3;

sigma_y1=12*my*h_sigma1/2./h_sigma1.^3;
sigma_y2=12*my*h_sigma2/2./h_sigma2.^3;

end