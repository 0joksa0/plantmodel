// #include <criterion/criterion.h>
// #include <solver.h>
// #include <stdio.h>
//
// #include "model/input.h"
// #include <math.h>
// #include <stddef.h> // offsetof
// #include <string.h>
//
// // === Tvoja postojeća infrastruktura ===
// extern void simulate_days(int days, Input* input);
// extern Input generate_input(void); // ako nemaš, koristi svoj init blok umesto ovog
//
// static inline real_t compute_rgr(real_t FW_start, real_t FW_end, real_t duration_days)
// {
//     return (RLOG(FW_end) - RLOG(FW_start)) / duration_days;
// }
//
// // === Fotoperiod scenariji (literaturne RGR vrednosti) ===
// typedef struct {
//     int hours;
//     real_t lit_rgr;
//     real_t lambda_c;
//     real_t lambda_sni;
//     real_t lambda_sb;
// } Scenario;
//
// static const Scenario SCENARIOS[] = {
//     { 4, REAL(0.0680), REAL(0.82), REAL(0.16), REAL(0.00344) },
//     { 6, REAL(0.1135), REAL(0.79), REAL(0.15), REAL(0.00413) },
//     { 8, REAL(0.1708), REAL(0.67), REAL(0.13), REAL(0.00515) },
//     { 12, REAL(0.2600), REAL(0.62), REAL(0.08), REAL(0.00578) },
//     { 18, REAL(0.3065), REAL(0.55), REAL(0.004), REAL(0.00524) },
// };
//
// static const int N_SCEN = sizeof(SCENARIOS) / sizeof(SCENARIOS[0]);
//
// // === Specifikacija polja koja menjamo (OAT) ===
// // Koristimo offsetof da bismo generički pristupili poljima iz Input.
// typedef struct {
//     const char* name;
//     size_t offset;
//     // opciono: ograničenja ili tip variranja (add vs mul). Ovde radimo multiplicativno.
// } FieldSpec;
//
// // Helper da dohvatimo pokazivač na polje (real_t*) po offsetu
// static inline real_t* field_ptr(Input* in, size_t off)
// {
//     return (real_t*)(((char*)in) + off);
// }
//
// // OVDE izlistaj polja koja želiš da “tweakuješ”.
// // Počni sa ključnim iz tvog modela; lako dodaš/oduzmeš.
// static const FieldSpec FIELDS[] = {
//     // {"lambda_sdr",                  offsetof(Input, lambda_sdr)},
//     // {"lambda_sdi",                  offsetof(Input, lambda_sdi)},
//     // {"lambda_sni",                  offsetof(Input, lambda_sni)},
//     { "lambda_sb", offsetof(Input, lambda_sb) },
// };
// static const int N_FIELDS = sizeof(FIELDS) / sizeof(FIELDS[0]);
//
// // Multiplikatori za OAT (možeš da proširiš/izmeniš)
// static const real_t MULTS[] = { REAL(1.21),REAL(1.23),REAL(1.25),REAL(1.27),REAL(1.29) };
// static const int N_MULTS = sizeof(MULTS) / sizeof(MULTS[0]);
//
// // Jedan run: vrati model RGR i apsolutnu grešku za dati fotoperiod i postavke
// static void run_one(const Scenario* sc, const Input* base, size_t field_off, real_t mult,
//     real_t* out_model_rgr, real_t* out_abs_err)
// {
//     Input in = *base; // kopija
//     // Primeni fotoperiod
//     in.photoperiod = REAL(sc->hours);
//     in.feedback_on_photosynthesis = sc->lambda_c;
//     in.lambda_sni = sc->lambda_sni;
//     in.lambda_sb = sc->lambda_sb;
//
//     // Primeni OAT promenu na izabrano polje
//     real_t* p = field_ptr(&in, field_off);
//     real_t original = *p;
//     *p = original * mult;
//
//     // Resetuj veličine koje treba da budu na početku simulacije (po potrebi)
//     // Ako imaš dnevni “trigger” ili noćni start škroba, setuj ovde.
//     real_t FW_start = in.total_biomass;
//
//     simulate_days(30, &in); // 30 dana
//
//     real_t model_rgr = compute_rgr(FW_start, in.total_biomass, REAL(30.0));
//     real_t abs_err = fabs(model_rgr - sc->lit_rgr);
//
//     if (out_model_rgr)
//         *out_model_rgr = model_rgr;
//     if (out_abs_err)
//         *out_abs_err = abs_err;
// }
//
// // Upis u CSV
// static FILE* open_csv(const char* path)
// {
//     FILE* f = fopen(path, "w");
//     if (!f) {
//         perror("fopen CSV");
//         return NULL;
//     }
//     fprintf(f, "field,multiplier,photoperiod_h,model_rgr,lit_rgr,abs_error\n");
//     return f;
// }
//
// // Glavna “sweep” rutina: za svako polje i multiplikator, prođi sve fotoperiode, zabeleži rezultate i agregiraj grešku
// static void param_sweep_to_csv(const char* csv_path)
// {
//     FILE* csv = open_csv(csv_path);
//     if (!csv)
//         return;
//
//     Input base = generate_input(); // ako nemaš, ovde ubaci tvoj “Input input = { ... }” blok
//
//     // Summary: čuvamo najbolji multiplikator po polju
//     real_t best_mult_per_field[N_FIELDS];
//     real_t best_mae_per_field[N_FIELDS];
//
//     for (int fi = 0; fi < N_FIELDS; ++fi) {
//         best_mult_per_field[fi] = MULTS[0];
//         best_mae_per_field[fi] = REAL(1e9);
//     }
//
//     for (int fi = 0; fi < N_FIELDS; ++fi) {
//         const FieldSpec* fs = &FIELDS[fi];
//
//         for (int mi = 0; mi < N_MULTS; ++mi) {
//             real_t mult = MULTS[mi];
//             real_t sum_abs_err = REAL(0.0);
//
//             for (int si = 0; si < N_SCEN; ++si) {
//                 const Scenario* sc = &SCENARIOS[si];
//                 real_t model_rgr, abs_err;
//                 run_one(sc, &base, fs->offset, mult, &model_rgr, &abs_err);
//
//                 sum_abs_err += abs_err;
//
//                 // CSV red
//                 fprintf(csv, "%s,%.6f,%d,%.6f,%.6f,%.6f\n",
//                     fs->name, (double)mult, sc->hours,
//                     (double)model_rgr, (double)sc->lit_rgr, (double)abs_err);
//             }
//
//             real_t mae = sum_abs_err / (real_t)N_SCEN;
//             if (mae < best_mae_per_field[fi]) {
//                 best_mae_per_field[fi] = mae;
//                 best_mult_per_field[fi] = mult;
//             }
//         }
//     }
//
//     fclose(csv);
//
//     // Štampaj sažetak najboljih multiplikatora
//     printf("\n=== OAT summary (best multiplier per field; metric: mean abs error across photoperiods) ===\n");
//     for (int fi = 0; fi < N_FIELDS; ++fi) {
//         printf("  %-30s  best_mult = %.3f   best_MAE = %.6f\n",
//             FIELDS[fi].name,
//             (double)best_mult_per_field[fi],
//             (double)best_mae_per_field[fi]);
//     }
//     printf("Results written to: %s\n", csv_path);
// }
//
// // Ako želiš kroz Criterion umesto main-a:
// Test(calibration, oat_param_sweep)
// {
//     param_sweep_to_csv("param_sweep.csv");
//     cr_assert(true); // nema faila; rezultat ide u CSV i stdout
// }
