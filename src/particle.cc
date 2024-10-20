#include "particle.h"

particle::particle(std::vector<double> &param, std::vector<int> &seq){
    for(int i=0;i<seq.size(); i++){
        if(seq[i] == 1) E0 = param[i];
        else if(seq[i] == 2) B0 = param[i] * pow(10., -9.);
        else if(seq[i] == 3) polarity = param[i];
        else if(seq[i] == 4) angle = param[i];
        else if(seq[i] == 5) D = param[i];
        else if(seq[i] == 6) indexA = param[i];
        else if(seq[i] == 7) Ek = param[i];
    }
}

particle::~particle(){}

const double particle::Wind(){
    double value;
    if(0.<theta && theta<pi/2.){
        value = 1.475 - 0.4 * tanh(6.8*((theta - 1.*pi/2.) + (15. + angle)/180.*pi)) * (3.5/5. - 1.5/5.*tanh((r-95.)/1.2));
    }
    else if(pi/2.<theta && theta<pi){
        value = 1.475 + 0.4 * tanh(6.8*((theta - 1.*pi/2.) - (15. + angle)/180.*pi)) * (3.5/5. - 1.5/5.*tanh((r-95.)/1.2));
    }

//     for(int i=0;i<9;i++){
//         double a = (5+10*i)*pi/180;
//         std::cout << a << "   " << 1.475 - 0.4 * tanh(6.8*((a - 1.*pi/2.) + (15. + angle)/180.*pi)) * (3.5/5. - 1.5/5.*tanh((r-95.)/1.2)) << std::endl;
//     }
//     for(int i=9;i<18;i++){
//         double a = (5+10*i)*pi/180;
//         std::cout << a << "   " << 1.475 + 0.4 * tanh(6.8*((a - 1.*pi/2.) - (15. + angle)/180.*pi)) * (3.5/5. - 1.5/5.*tanh((r-95.)/1.2)) << std::endl;
//     }
// getchar();
    return value * 400.;
}

const double particle::Theta_S(){
    double value;
    value = asin(sin(angle*pi/180.)*sin(Omega*r*AU/400.)/0.8354);
// std::cout << "theta_s :  " << value << "   " << sin(Omega*r*AU/400.) << std::endl;
// getchar();
    return value;
}

const double particle::Heav(){
    double theta_s = Theta_S();
    double value;
    if(theta<pi/2.-theta_s) value = 1.;
    else if(pi/2.-theta_s<theta) value = -1.;

    return value;
}

const double particle::B_r(const double &heaviside){
    double value;
    value = B0 * heaviside * polarity / pow(r, 2.);

    return value;
}

const double particle::B_p(const double &heaviside){
    double value;
    value = -1. * B0 * r * Omega *sin(theta) * heaviside * polarity / Vs;

    return value;
}

const double particle::K_rr(){
    double kx = k_n * D * pow(Ek, indexA);
    double ky = 0.02 * kx;
    double kr = kx * pow(cos(psy), 2.) + ky * pow(sin(psy), 2.);

    //     std::cout << "VD :  " << kx << "  " << D << "  " << ky << "  " <<  pow(Ek, indexA)<< std::endl;
    // getchar();

    return kr;
}

const double particle::K_tt(){
    double kx = k_n * D * pow(Ek, indexA);
    double ky = 0.02 * kx;
    double kz;
    if(theta < pi/2. && 0 < theta){
        kz = ky * (2. - 1.*tanh(8*((theta + (- 90. + 35.)*pi/180.))));
    }
    else if(pi/2.<theta && theta<pi){
        kz = ky * (2. + 1.*tanh(8*((theta + (- 90. - 35.)*pi/180.))));
    }

// std::cout << "VD :  " << ky<< "  " << kz << "  " << (2. + 1.*tanh(8*((theta + (- 90. - 35.)*pi/180.)))) << std::endl;


//         getchar();
    return kz;
}

const double particle::K_pp(){
    double kx = k_n * D * pow(Ek, indexA);
    double ky = 0.02 * kx;

    return ky;
}



void particle::step(){
    std::random_device rd;
    std::mt19937 gen(rd());
    double mean = 0.0;
    double dev = 1.0;
    std::normal_distribution<double> dist(mean, dev);

    double record_T = 0.;

    while(r<boundary || record_T<pow(10.,10.)){
        record_T += dt;

        rigidity = A/Z * sqrt(Ek*(Ek + 2.*E0));
        M_p = sqrt(Ek*(Ek + 2.*E0));
        V_p = sqrt(1. - 1./pow(Ek/E0+1., 2.)) * light;
        Vs = Wind()/AU;

        heaviside = Heav();
        Br = B_r(heaviside);
        Bp = B_p(heaviside);
        psy = atan(fabs(Bp/Br)); 
        
        
        k_rr = K_rr();
        k_tt = K_tt();
        k_pp = K_pp();


        double gamma = r*Omega*sin(theta)/Vs;
        double drift = 2.*M_p*pow(10.,6.)*V_p/light*r*AU/3./(4.7/5.*pow(10.,21.)*B0)/AU;
        Vdr_gc = drift/pow(1+gamma*gamma, 2.) * heaviside * (-1.*gamma/ tan(theta));
        Vdp_gc = drift/pow(1+gamma*gamma, 2.) * heaviside * (2. + gamma*gamma) * gamma;
        Vdt_gc = drift/pow(1+gamma*gamma, 2.) * heaviside * gamma*gamma / tan(theta);
// std::cout << "VD :  " << drift<< "  " << Vdr_gc << "  " << Vdp_gc << "  " << Vdt_gc << std::endl;
//         getchar();
        double Vd = drift/(1+gamma*gamma) ;




        double L0, Rg;


        double dw = gen();

        r += (-1.*Vs - Vdr_gc - Vdr_HCS + 2./r) * dt
            + sqrt(2. * k_rr * dt) * dw;

        theta += (-1.*Vdt_gc/r + 1./(r*r*sin(theta))*cos(theta)*k_tt) * dt 
                + 1./r * sqrt(2.*k_tt*dt) * dw;

        phy += (-1.*Vdp_gc - Vdp_HCS) / (r * sin(theta)) * dt
                + sqrt(2.*k_pp*dt) * dw / (r * sin(theta));

        Ek += 2.*V_p / (3.*r) * (Ek*Ek + 2.*Ek*E0) / (Ek + E0) * dt;

        if(r<0.) {r = 0.; break;}

        if(theta<0.) {theta = fabs(theta); phy += pi;}
        else if(pi<theta) {theta = 2.*pi - theta; phy += pi;}

        if(phy<0.) phy = 2.*pi - phy;
        else if(2.*pi<phy) phy -= 2.*pi;



    }

}