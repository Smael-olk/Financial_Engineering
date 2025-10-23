#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>



struct OptionData {
    double S0;    // spot price
    double K;     // strike
    double r;     // risk-free rate
    double d;     // dividend yield
    double vol;   // volatility
    double T;     // time to maturity
};


class OptionPricer{
    public :
        virtual double price() const =0;
        virtual ~OptionPricer() = default;
};

class blsPricer : public OptionPricer{
    private :
        OptionData opt;
    public : 
        explicit blsPricer(const OptionData& data): opt(data){}
        double price() const override{
            double d1=(std::log(opt.S0/opt.K)+(opt.r-opt.d+0.5*opt.vol*opt.vol)*opt.T)/(opt.vol*std::sqrt(opt.T));
            double d2=d1-opt.vol*std::sqrt(opt.T);
            double Nd1=0.5*(1.0+std::erf(d1/std::sqrt(2.0)));
            double Nd2=0.5*(1.0+std::erf(d2/std::sqrt(2.0)));
            return opt.S0*std::exp(-opt.d*opt.T)*Nd1-opt.K*std::exp(-opt.r*opt.T)*Nd2;
        }
};

class crrPricer : public OptionPricer {
private:
    OptionData opt;
    int steps;
public:
    crrPricer(const OptionData& data, int n) : opt(data), steps(n) {}

    double price() const override {
        double dt = opt.T / steps;
        double u = std::exp(opt.vol * std::sqrt(dt));
        double d = 1.0 / u;
        double p = (std::exp((opt.r - opt.d) * dt) - d) / (u - d);

        std::vector<double> prices(steps + 1);
        for (int i = 0; i <= steps; ++i) {
            double ST = opt.S0 * std::pow(u, steps - i) * std::pow(d, i);
            prices[i] = std::max(ST - opt.K, 0.0); 
        }


        for (int j = steps - 1; j >= 0; --j) {
            for (int i = 0; i <= j; ++i) {
                prices[i] = std::exp(-opt.r * dt) * (p * prices[i] + (1 - p) * prices[i + 1]);
            }
        }

        return prices[0];
    }
};


class mcPricer  : public OptionPricer {
    private : 
        OptionData opt;
        int nPaths;
    public :
        mcPricer(const OptionData& data, int paths ) : opt(data),nPaths(paths){}
        double price() const override{
            std::mt19937_64 rng(42);
            std::normal_distribution<double> norm(0.0,1.0);
            double sum_payoffs=0.0;
            for(int i=0;i<nPaths;i++){
                double Z=norm(rng);
                double ST=opt.S0*std::exp((opt.r-opt.d-0.5*opt.vol*opt.vol)*opt.T+opt.vol*std::sqrt(opt.T)*Z);
                double payoff=std::max(0.0,ST-opt.K);
                sum_payoffs+=payoff;
            }
            double discounted_payoff=std::exp(-opt.r*opt.T)*(sum_payoffs/nPaths);
            return discounted_payoff;
        }
};

int main() {
    OptionData data_FE{1.00, 1.05, 0.025, 0.02, 0.21, 1.0 / 3.0};
    OptionData data_MSTR{1.00, 1.05, 0.025, 0.02, 0.21, 1.0 / 3.0};

    blsPricer bs(data_FE);
    crrPricer crr(data_FE, 50);
    mcPricer mc(data_FE, 100000);

    std::cout << "Black Scholes price: " << bs.price() << "\n";
    std::cout << "CRR Tree price:      " << crr.price() << "\n";
    std::cout << "Monte Carlo price:   " << mc.price() << "\n";

    return 0;
}
