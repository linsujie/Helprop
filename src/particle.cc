#include "particle.h"

using namespace std;
using namespace Unit;

particle::particle() {}

particle::particle(const map<string, docopt::value> &args){
    auto fargs = [&](const string& key) -> double { return stod(args.at(key).asString()); };
    mass = fargs("--mass") * GeV;
    B0 = fargs("--B0") * nT;
    polarity = args.at("--polarity").asLong();
    angle = fargs("--angle") * deg;
    D = fargs("--D") * 1e22 * cm * cm / sec;
    indexA = fargs("--indexA");

    Ek = mass;
}

particle::~particle(){}

const double particle::Wind(){
    double value;
    if(0.<theta && theta<pi/2.){
        value = 1.475 - 0.4 * tanh(6.8*((theta - pi/2) + (15*deg + angle))) * (3.5/5. - 1.5/5.*tanh((r-95*AU)/1.2/AU));
    }
    else if(pi/2.<theta && theta<pi){
        value = 1.475 + 0.4 * tanh(6.8*((theta - pi/2) - (15*deg + angle))) * (3.5/5. - 1.5/5.*tanh((r-95*AU)/1.2/AU));
    }

//     for(int i=0;i<9;i++){
//         double a = (5+10*i)*pi/180;
//         cout << a << "   " << 1.475 - 0.4 * tanh(6.8*((a - 1.*pi/2.) + (15. + angle)/180.*pi)) * (3.5/5. - 1.5/5.*tanh((r-95.)/1.2)) << endl;
//     }
//     for(int i=9;i<18;i++){
//         double a = (5+10*i)*pi/180;
//         cout << a << "   " << 1.475 + 0.4 * tanh(6.8*((a - 1.*pi/2.) - (15. + angle)/180.*pi)) * (3.5/5. - 1.5/5.*tanh((r-95.)/1.2)) << endl;
//     }
// getchar();
    return value * 400 * (km/sec);
}

const double particle::Theta_S(){
    double value;
    value = asin(sin(angle)*sin(Omega*r/(400*km/sec))/0.8354);
// cout << "theta_s :  " << value << "   " << sin(Omega*r/(400*km/sec))) << endl;
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
    const double r0 = 1 * AU;
    return B0 * heaviside * polarity / pow(r0/r, 2.);
}

const double particle::B_p(const double &heaviside){
    return -1. * B0 * r * Omega *sin(theta) * heaviside * polarity / Vs;
}

const double particle::K_rr(){
    double kx = D * pow(Ek/GeV, indexA);
    double ky = 0.02 * kx;
    double kr = kx * pow(cos(psi), 2.) + ky * pow(sin(psi), 2.);

    //     cout << "VD :  " << kx << "  " << D << "  " << ky << "  " <<  pow(Ek, indexA)<< endl;
    // getchar();

    return kr;
}

const double particle::K_tt(){
    double kx = D * pow(Ek/GeV, indexA);
    double ky = 0.02 * kx;
    double kz;
    if(theta < pi/2. && 0 < theta)
        kz = ky * (2. - 1.*tanh(8*((theta + (- 90. + 35.)*deg))));
    else if(pi/2.<theta && theta<pi)
        kz = ky * (2. + 1.*tanh(8*((theta + (- 90. - 35.)*deg))));

// cout << "VD :  " << ky<< "  " << kz << "  " << (2. + 1.*tanh(8*((theta + (- 90. - 35.)*pi/180.)))) << endl;


//         getchar();
    return kz;
}

const double particle::K_pp(){
    double kx = D * pow(Ek/GeV, indexA);
    double ky = 0.02 * kx;

    return ky;
}

double particle::get_HCS_distance() const {
    return 0;
}

void particle::step() {
    random_device rd;
    mt19937 gen(rd());
    double mean = 0.0;
    double dev = 1.0;
    normal_distribution<double> dist(mean, dev);

    double record_T = 0.;

    while(r<boundary || record_T<pow(10.,10.)){
        record_T += dt;

        M_p = sqrt(Ek*(Ek + 2*mass));
        rigidity = A/Z * M_p;
        V_p = M_p/(Ek + mass) * light;
        Vs = Wind();

        heaviside = Heav();
        Br = B_r(heaviside);
        Bp = B_p(heaviside);
        psi = atan(fabs(Bp/Br)); 
        
        
        k_rr = K_rr();
        k_tt = K_tt();
        k_pp = K_pp();

        double r0 = 1.0 * AU;
        double gamma = r*Omega*sin(theta)/Vs;
        double drift = 2 * M_p * V_p * r / (3 * Z * e * B0 * r0 * r0 * light);
         //2.*(M_p/MeV)*pow(10.,6.)*(V_p/(AU/sec))/(light/(AU/sec))*(r/AU)/3./(4.7/5.*pow(10.,21.)*B0);
         // 我用统一单位计算之后，和你的表达结果对不上，估计我们要讨论一下

        Vdr_gc = drift/pow(1+gamma*gamma, 2.) * heaviside * (-1.*gamma/ tan(theta));
        Vdp_gc = drift/pow(1+gamma*gamma, 2.) * heaviside * (2. + gamma*gamma) * gamma;
        Vdt_gc = drift/pow(1+gamma*gamma, 2.) * heaviside * gamma*gamma / tan(theta);
        //cout << "Z e B0: " << Z << " " <<  e << " " << B0 << " " << r / AU << " " << r0 / AU << " " << light << " " << k_rr / (cm * cm / sec) << " " << (-1.*Vs - Vdr_gc - Vdr_HCS) * dt / AU << " " << sqrt(k_rr * dt) / (AU) << endl;
        //cout << "momentum : " << M_p / GeV << " " << V_p / (km/sec) << " " << Z << " " << Z << " " << A << endl;
        //cout << "VD :  " << drift / (km/sec) << "  " << Vdr_gc / (km/sec) << "  " << Vdp_gc / (km/sec) << "  " << Vdt_gc / (km/sec) << endl;
//         getchar();
        double Vd = drift/(1+gamma*gamma) ;



        double L0, Rg;


        double dw = dist(gen);
        double d_HCS = get_HCS_distance(); // 这个函数我之后填，目前我连HCS的函数形式都没找着

        //r += (-1.*Vs - Vdr_gc - Vdr_HCS + 2./r) * dt
        cout << "rbefore: " << r / AU << endl;
        r += (-1.*Vs - Vdr_gc - Vdr_HCS) * dt
            + sqrt(2. * k_rr * dt) * dw;
        cout << "rafter: " << r / AU << " " << (-1.*Vs - Vdr_gc - Vdr_HCS) * dt / AU << " " << sqrt(2. * k_rr * dt) / AU << " " << dw << endl;

        theta += (-1.*Vdt_gc/r + 1./(r*r*sin(theta))*cos(theta)*k_tt) * dt 
                + 1./r * sqrt(2.*k_tt*dt) * dw;

        phi += (-1.*Vdp_gc - Vdp_HCS) / (r * sin(theta)) * dt
                + sqrt(2.*k_pp*dt) * dw / (r * sin(theta));

        Ek += 2.*V_p / (3.*r) * (Ek*Ek + 2.*Ek*mass) / (Ek + mass) * dt;

        if(r<0.) {r = 0.; break;}

        if(theta<0.) {theta = fabs(theta); phi += pi;}
        else if(pi<theta) {theta = 2.*pi - theta; phi += pi;}

        if(phi<0.) phi = 2.*pi - phi;
        else if(2.*pi<phi) phi -= 2.*pi;
    }

}