#include <schrodinger.h>

int main(int argc, char **argv) {
    /*
        Choix des valeurs => ordre de grandeur proche de 1.
        Unités des grandeurs physiques : 
        m       1e-30 kg
        L       nm (1e-9 m)
        E       eV (1e-19 J)
        hbar    1.05e-34 J.s 

        Pour simplifier les ordres de grandeurs, on prend :

        hbarc   en eV*nm, au lieu de hbar (c : célérité de la lumière)
        hbarc = 197 eV*nm;

        m       en MeV/c^2, au lieu de kg
        m_electron = 9.11e-31 kg = 0.511 MeV/c^2 = 511000 eV/c^2
    */

    double m = 511000.0;   //Masse de l'electron, en eV/c^2
    double L = 1.0;       //Largeur du puit, en nm (1e-9 m)
    double E = 0.0;       //Energie de la particule, en eV

    // Resolution : Cas du puit de potentiel infini
    puit_potentiel_infini(m, L, E);
    return 0;
}