# Data manifest

Every data file read by the code, relative to the data directory. Paths match the code, so files from the CSL site or other sources should be placed at the same relative path.

Status:
- **in repo**: committed under `data/` (SHA-256 given).
- **CSL site**: too large for GitHub; available from https://csl.noaa.gov/groups/csl8/modeldata/.
- **external**: third-party data not redistributed here; obtain from the listed source.
- **not in repo**: needed only by the original programs in `code/original/`; available on request.

"Used by" says whether the streamlined programs in `code/steps/` read the file, or only the original programs in `code/original/`. The steps write their own intermediate files as NetCDF in `output/intermediate/`, so the "intermediate" save files listed here are needed only by the original programs.

| Path | Status | Size (MB) | Type | Used by | Source | SHA-256 |
|---|---|---|---|---|---|---|
| `Aircraft/Airborne_Save_Files/Aircraft_mean_ages_early_bins_comb.sav` | in repo | 0.0 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `59fbad254d0ab5fb65da89032696805ce5cb4217d9124e95fc7cb6dd318ffe9d` |
| `Aircraft/Airborne_Save_Files/Andrews_N2O_age.sav` | external | 0.0 | input | code/steps | Andrews et al. N2O vs. mean age relationship; available from the authors |  |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_early_bins.sav` | CSL site | 53.7 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_late_bins.sav` | in repo | 40.4 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `1e549dc7e93dfdab13886c25da532eb5b892854de1a9c67306eb362e5efd3270` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_mid.sav` | in repo | 43.4 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `2471415f0cead96fbbac157d486e313036d3233a246576db0633445781fbeb1b` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_mid_bins.sav` | in repo | 46.3 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `e3686392bef86f8bc661070082e02700acae035c861f70b68e087c54149276cc` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_optimum_sf6_co2_mid_bins.sav` | CSL site | 3927.2 | input | code/original only | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_early_bins.sav` | in repo | 22.6 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `d36cf5cb3807897ec8219480b7073683dd50e487b02546822414c9615a490b04` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_late.sav` | in repo | 0.9 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `5a70310e8ea6fd6359aeadced90e21f552b119c686ccb469f0f32a9ac08db93e` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_late_bins.sav` | in repo | 22.6 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `fd9f120989b2e728d3c61e5a8f8fade68aaa479d01c92629196de75906b59896` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_mid.sav` | in repo | 4.0 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `d59079c6a652dc75f3ab20be09b1c29d9e4361fc1daf9259c7fa9361dfbcc06f` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_mid_bins.sav` | in repo | 23.2 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `5c83dc0cb7d03acb38cd38b62e0e16c47291b7d227618dcaf468a9ed49b95e04` |
| `Aircraft/Airborne_Save_Files/aircraft_mean_ages_sweep_sf6_co2_early.sav` | CSL site | 344.9 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Airborne_Save_Files/ashoe.sav` | CSL site | 2293.6 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/Mean_age_relationships_aircraft.sav` | not in repo | 0.1 | input | code/original only | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/N2O_CH4_norm_curve_post2000.sav` | not in repo | 0.0 | input | code/original only | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_age_grid_early.sav` | CSL site | 65.9 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_age_grid_late.sav` | CSL site | 65.3 | input | code/original only | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_early.sav` | CSL site | 78.3 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_late.sav` | in repo | 15.5 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `9d2d6583c82a1788c2ee7d300708bdcb22af0acafadab86f1b758fcade225004` |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_early.sav` | CSL site | 7988.0 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_late.sav` | CSL site | 1770.5 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_mid.sav` | CSL site | 4576.4 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_optimum_sf6_co2_early_bins.sav` | CSL site | 3927.2 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_optimum_sf6_co2_late_bins.sav` | CSL site | 3927.2 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_sf6_early.sav` | in repo | 0.7 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) | `26f950d2f36b55aeab6f2f59651d59134afffcd6ca2e0a35d49a0212e4457864` |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_sweep_sf6_co2_late.sav` | CSL site | 317.7 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Aircraft/Missions/idlsave_files/aircraft_mean_ages_sweep_sf6_co2_mid.sav` | CSL site | 53.5 | input | code/steps | Merged from aircraft mission archives (ESPO, DCOTSS, SABRE; see README) |  |
| `Balloon/Age_N2O_gridded_N2O_elat_adj_time_series.sav` | CSL site | 1512.2 | intermediate | code/original only | This study (merged / derived) |  |
| `Balloon/Aircore_common_merge.sav` | in repo | 10.4 | input | code/steps | Merged from NOAA GML AirCore, https://gml.noaa.gov/ccgg/arc/?id=144 | `41c30ce0ee8249b6e33d8dca6dab854ac5c66947f7980e6df92da3a415f63fba` |
| `Balloon/Aircore_grid_elat_adj.sav` | not in repo | 9.6 | intermediate | code/original only | Merged from NOAA GML AirCore, https://gml.noaa.gov/ccgg/arc/?id=144 |  |
| `Balloon/Aircore_n2o_ch4_v2024.sav` | in repo | 0.0 | input | code/steps | Merged from NOAA GML AirCore, https://gml.noaa.gov/ccgg/arc/?id=144 | `9ba81f0dcf37d84de3ffc597033d280fc2ff690ace79f2a2c7760908c0bf692e` |
| `Balloon/Balloon_Tracers_vs_mean_age.sav` | not in repo | 0.0 | input | code/original only | This study (merged / derived) |  |
| `Balloon/Balloon_mean_ages_co2_early.sav` | in repo | 3.2 | input | code/steps | This study (merged / derived) | `2b6b7090d640e4aa0dd6139fc9e299a6c8d98021f9246d13a510df6d1d8521b6` |
| `Balloon/Balloon_mean_ages_co2_early_bins.sav` | in repo | 37.8 | input | code/steps | This study (merged / derived) | `3b5de9ead271d512dba2f91d4a29ec09bad51a49fa7bd9042b3d550a602351b8` |
| `Balloon/Balloon_mean_ages_co2_late_v2025.sav` | in repo | 7.9 | input | code/steps | This study (merged / derived) | `e3c366bc405151eb9fdd571b42b7debdacbcf94d8c43fc738126fcdf9cca3b44` |
| `Balloon/Balloon_mean_ages_co2_late_v2025_bins.sav` | in repo | 39.7 | input | code/steps | This study (merged / derived) | `3a6732b97804388b24db44fefccb4b301b54e12cd3b18935b02abdd9b15223dc` |
| `Balloon/Balloon_mean_ages_optimum_sf6_co2_early_bins.sav` | CSL site | 645.3 | input | code/original only | This study (merged / derived) |  |
| `Balloon/Balloon_mean_ages_sf6_early.sav` | in repo | 0.2 | input | code/steps | This study (merged / derived) | `03f35688aa984a84fc81f2029eb2629a9b429f338ff6d27a35e3b44d16c7ebc7` |
| `Balloon/Balloon_mean_ages_sf6_early_bins.sav` | in repo | 22.3 | input | code/steps | This study (merged / derived) | `182b409f1349477626d6f5b6c88d3b57d3ddb2fc0f61b42761b437bf51679760` |
| `Balloon/Balloon_mean_ages_sf6_late.sav` | in repo | 0.1 | input | code/steps | This study (merged / derived) | `6ec2ea4ecb1be2838772da9223679898d402ed2addd40ec045653a83bd4d3a5e` |
| `Balloon/Balloon_mean_ages_sf6_late_bins.sav` | in repo | 24.8 | input | code/steps | This study (merged / derived) | `11278df05f17899b99eb636be9dacfa3ffc8f45c69b47e0c4a03620596ef48f2` |
| `Balloon/Balloon_mean_ages_sf6_late_v2025_UTC.sav` | not in repo | 0.0 | input | code/original only | This study (merged / derived) |  |
| `Balloon/Balloon_mean_ages_sweep_sf6_co2_early.sav` | in repo | 34.5 | input | code/steps | This study (merged / derived) | `cb4f44c9c62823b8e18ca30a9fd56198a1642a088dd4c8cd96587ffc4ff875f2` |
| `Balloon/Balloon_mean_ages_sweep_sf6_co2_late_v2025.sav` | in repo | 31.1 | input | code/steps | This study (merged / derived) | `3952c3886e52453b4d60a699a4d00991b633e85106126f5a5a0db91dd98932ae` |
| `Balloon/Elat_adj_all.sav` | not in repo | 39.4 | intermediate | code/original only | This study (merged / derived) |  |
| `Balloon/Engel_balloon_all5.txt` | external | 2.1 | input | code/original only | Compiled cryo-flask profiles from Engel et al. (2009), Ray et al. (2014) and Sugawara et al. (2025); contact the data providers |  |
| `Balloon/Engel_mean_ages_opt.sav` | external | 6.8 | input | code/original only | Cryo-flask balloon data from Engel et al. (paper refs 9, 52); contact the data providers |  |
| `Balloon/Engel_profiles_v2025.sav` | external | 22.7 | intermediate | code/steps | Combined flask, OMS, AirCore and WAS profiles built by age_time_series.pro (read1-read5 = 1) from Engel et al. data; contact the data providers |  |
| `Balloon/Mean_age_and_N2O_trend_alt_profile_elat_adj.sav` | not in repo | 17.1 | intermediate | code/original only | This study (merged / derived) |  |
| `Balloon/Mean_age_relationships.sav` | in repo | 0.1 | input | code/steps | This study (merged / derived) | `59c7f2b206803837f1ca71babd5b1f3fe92b5a3dc2f49fb7a06f544ac5dc7708` |
| `Balloon/Mean_age_vs_N2O_1990s.sav` | not in repo | 0.1 | intermediate | code/original only | This study (merged / derived) |  |
| `Balloon/NSamp_OMS.sav` | in repo | 0.0 | intermediate | code/steps | Merged from OMS balloon data, https://espoarchive.nasa.gov/archive/browse | `daaa43f5cd748a89030d5e86003e805764912fb59129cd13ac6ffc5eb567ad76` |
| `Balloon/OMS_common_gc.sav` | not in repo | 46.7 | input | code/original only | Merged from OMS balloon data, https://espoarchive.nasa.gov/archive/browse |  |
| `Balloon/OMS_common_gc_merge.sav` | in repo | 19.0 | input | code/steps | Merged from OMS balloon data, https://espoarchive.nasa.gov/archive/browse | `7465b001c161aff9dbe881dad2119b374cafc989ca95336c9522bcdcf86c24b8` |
| `Balloon/OMS_grid_elat_adj.sav` | not in repo | 10.6 | intermediate | code/original only | Merged from OMS balloon data, https://espoarchive.nasa.gov/archive/browse |  |
| `Balloon/OMS_tracer_profile_avgs.sav` | in repo | 0.0 | input | code/steps | Merged from OMS balloon data, https://espoarchive.nasa.gov/archive/browse | `6bfc2ce5cd29c4fe2ff3d055f1dcb3cf2e2618e73f09379fc8fca1194d8ad137` |
| `Balloon/Sample_type_time_series.sav` | not in repo | 0.0 | intermediate | code/original only | This study (merged / derived) |  |
| `Balloon/WAS_mean_ages.sav` | not in repo | 0.0 | input | code/original only | This study (merged / derived) |  |
| `Balloon/balloon_mean_age_grid_early.sav` | CSL site | 58.9 | input | code/original only | This study (merged / derived) |  |
| `Balloon/balloon_mean_age_grid_late_v2025.sav` | CSL site | 160.8 | input | code/original only | This study (merged / derived) |  |
| `Balloon/balloon_mean_age_seas_grid_early.sav` | not in repo | 0.1 | intermediate | code/original only | This study (merged / derived) |  |
| `Models/CCMI/CCMI-1_refc1_means.sav` | external | 552.7 | input | code/original only | CCMI-2022 / CCMI-1 output via CEDA, https://data.ceda.ac.uk/badc/ccmi/data/post-cmip6/ccmi-2022 |  |
| `Models/CCMI/CCMI-2022_refd1_means.sav` | external | 2927.7 | input | code/steps | CCMI-2022 / CCMI-1 output via CEDA, https://data.ceda.ac.uk/badc/ccmi/data/post-cmip6/ccmi-2022 |  |
| `Models/CCMI/CCMI-2022_refd1_trends.sav` | external | 7.9 | input | code/steps | CCMI-2022 / CCMI-1 output via CEDA, https://data.ceda.ac.uk/badc/ccmi/data/post-cmip6/ccmi-2022 |  |
| `Models/CCMI/CCMI-2022_refd2_means.sav` | not yet available |  | input | code/original only | CCMI-2022 output via CEDA, https://data.ceda.ac.uk/badc/ccmi/data/post-cmip6/ccmi-2022 |  |
| `Models/CCMI/CCMI-2022_refd2_trends.sav` | external | 15.6 | input | code/steps | CCMI-2022 / CCMI-1 output via CEDA, https://data.ceda.ac.uk/badc/ccmi/data/post-cmip6/ccmi-2022 |  |
| `Models/CLaMS/CLaMS_ERA5_N2O_on_mean_ages.sav` | external | 17.0 | input | code/original only | CLaMS model output; available from the model developers |  |
| `Models/CLaMS/CLaMS_SD_UVT_ERA5_tracers.sav` | external | 113.8 | input | code/original only | CLaMS model output; available from the model developers |  |
| `Models/EMAC/SF6_age_bias_correction_Garny.txt` | external | 0.0 | input | code/original only | SF6 age bias correction (Garny); available from the authors |  |
| `Models/TLP/ideal_sweep_w_LowStrat.sav` | in repo | 13.5 | input | code/steps | Tropical Leaky Pipe model output (this study) | `f5532d126b92c372196d66bed48322b40e0725066e500a6ec844ec7b91cb7c10` |
| `Models/TLP/ideal_sweep_w_LowStrat_UpStrat_Budget.sav` | CSL site | 276.8 | input | code/steps | Tropical Leaky Pipe model output (this study) |  |
| `Models/TLP/ideal_sweep_w_UpStrat_budget.sav` | CSL site | 66.5 | input | code/steps | Tropical Leaky Pipe model output (this study) |  |
| `Models/WACCM/FWSD_early_aircraft_balloon_subsample_mean_age_N2O.sav` | external | 0.0 | input | code/original only | WACCM FWSD output; available from the model developers |  |
| `Models/WACCM/FWSD_means.sav` | external | 296.3 | input | code/steps | WACCM FWSD output; available from the model developers |  |
| `Models/WACCM/WACCM_FWSD_N2O_on_age.sav` | external | 3.4 | input | code/original only | WACCM FWSD output; available from the model developers |  |
| `Reanalysis/MERRA2/tp.monmean.zm.nc` | in repo | 5.1 | input | code/steps | Derived from NASA MERRA-2, https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/ | `3c21657c703dc09f5d3911e542815fdd41c42e15bdc8d24ec2f134e51dffc230` |
| `Satellite/ACE/ACE_gridding_5p3.sav` | external | 879.0 | input | code/steps | ACE-FTS (registration required), https://databace.scisat.ca/ |  |
| `Surface_Trace_Gas/CH4/CH4_mbl_30S-30N.txt` | not in repo | 0.1 | input | code/original only | NOAA GML Marine Boundary Layer Reference, https://gml.noaa.gov/ccgg/mbl/ |  |
| `Surface_Trace_Gas/CO2/CO2_mlo_monthly.txt` | in repo | 0.0 | input | code/steps | NOAA GML Mauna Loa / Samoa CO2, https://gml.noaa.gov/ccgg/trends/ | `a9ef9ab9bf1270e5ee38c3ec1f81d5ba4dfd51ee87d59f5363216f53608db6ad` |
| `Surface_Trace_Gas/CO2/CO2_smo_monthly.txt` | in repo | 0.0 | input | code/steps | NOAA GML Mauna Loa / Samoa CO2, https://gml.noaa.gov/ccgg/trends/ | `f40d9d8cfbb794bae487e3436e69de93db2bad642ce32276b9b9544077e7d4ac` |
| `Surface_Trace_Gas/GMD/cfc12/combined/HATS_global_F12.txt` | in repo | 0.2 | input | code/steps | NOAA GML HATS combined CFC-12, https://gml.noaa.gov/hats/combined/CFC12.html | `1a5c51d28b9f924a1d1fa0ad0cb07d571a37f88836212bf91125838949f99993` |
| `Surface_Trace_Gas/GMD/n2o/combined/GML_global_N2O.txt` | in repo | 0.2 | input | code/steps | NOAA GML HATS combined N2O, https://gml.noaa.gov/hats/combined/N2O.html | `5094bc7f7eab00360f7c3bf54059b422d9e701250b634a1ff669b202a2e4d368` |
| `Surface_Trace_Gas/Tracer_entry_time_series.sav` | CSL site | 582.7 | input | code/original only | Derived from NOAA GML surface data |  |
| `swoosh/n2o_merge.sav` | in repo | 12.0 | input | code/steps | Derived from NOAA SWOOSH, https://csl.noaa.gov/groups/csl8/swoosh/ | `03a69ace3f6bb60d3d65fbac6b26735defa4a21e582d3980a3a3dd44174d3b31` |
