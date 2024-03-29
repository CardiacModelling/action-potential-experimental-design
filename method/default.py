###############################################################################
## Defining Model Conversion, conductance, ion conc etc.
###############################################################################

# 'ina.s', 'ical.s' 'ikr.s', 'iks.s', 'ito.s', 'inaca.s', 'ik1.s'
# 'inak.s', 'if.s'

# model_voltage = {'mmt_file_name': voltage_variable_name_in_mmt}
model_voltage = {
    'tnnp-2004.mmt': 'membrane.V',
    'fink-2008.mmt': 'membrane.V',
    'fink-2008-lei-2019-herg.mmt': 'membrane.V',
    'fink-2008-lei-2019-herg-tunable.mmt': 'membrane.V',
    'grandi-2010.mmt': 'membrane_potential.V_m',
    'ohara-2011.mmt': 'membrane.V',
    'cipa-2017.mmt': 'membrane.v',
    'tomek-2019.mmt': 'membrane.v',
}

model_pace = {
    'tnnp-2004.mmt': 'engine.pace',
    'fink-2008.mmt': 'engine.pace',
    'fink-2008-lei-2019-herg.mmt': 'engine.pace',
    'fink-2008-lei-2019-herg-tunable.mmt': 'engine.pace',
    'grandi-2010.mmt': 'interface.pace',
    'ohara-2011.mmt': 'engine.pace',
    'cipa-2017.mmt': 'membrane.pace',
    'tomek-2019.mmt': 'membrane.pace',
}

# model_current = {'mmt_file_name': current_dictionary}
# current_dictionary = { conductance_scaling: current_in_mmt }
# different model might call current differently...
# this was mainly for voltage clamp simulation...

tnnp_2004_current = {
    'ina.s': 'fast_sodium_current.i_Na',
    'ical.s': 'L_type_Ca_current.i_CaL',
    'ikr.s': 'rapid_time_dependent_potassium_current.i_Kr',
    'iks.s': 'slow_time_dependent_potassium_current.i_Ks',
    'ito.s': 'transient_outward_current.i_to',
    'inaca.s': 'sodium_calcium_exchanger_current.i_NaCa',
    'ik1.s': 'inward_rectifier_potassium_current.i_K1',
    'inak.s': 'sodium_potassium_pump_current.i_NaK',
}

fink_2008_current = {
    'ina.s': 'ina.i_Na',
    'ical.s': 'ical.i_CaL',
    'ikr.s': 'ikr.i_Kr',
    'iks.s': 'iks.i_Ks',
    'ito.s': 'ito.i_to',
    'inaca.s': 'inaca.i_NaCa',
    'ik1.s': 'ik1.i_K1',
    'inak.s': 'inak.i_NaK',
}

grandi_2010_current = {
    'ina.s': 'I_Na.I_Na',
    'ical.s': 'I_Ca.I_Ca',
    'ikr.s': 'I_Kr.I_kr',
    'iks.s': 'I_Ks.I_ks',
    'ito.s': 'I_to.I_to',
    'inaca.s': 'I_NCX.I_ncx',
    'ik1.s': 'I_Ki.I_ki',
    'inak.s': 'I_NaK.I_nak',
}

ohara_2011_current = {
    'ina.s': 'ina.INa',
    'ical.s': 'ical.ICaL',
    'ikr.s': 'ikr.IKr',
    'iks.s': 'iks.IKs',
    'ito.s': 'ito.Ito',
    'inaca.s': 'inaca.INaCa',
    'ik1.s': 'ik1.IK1',
    'inak.s': 'inak.INaK',
}

cipa_2017_current = {
    'ina.s': 'INa.INa',
    'ical.s': 'ICaL.ICaL',
    'ikr.s': 'IKr.IKr',
    'iks.s': 'IKs.IKs',
    'ito.s': 'Ito.Ito',
    'inaca.s': 'INaCa_i.INaCa',
    'ik1.s': 'IK1.IK1',
    'inak.s': 'INaK.INaK',
}

tomek_2019_current = {
    'ina.s': 'INa.INa',
    'ical.s': 'ICaL.ICaL',
    'ikr.s': 'IKr.IKr',
    'iks.s': 'IKs.IKs',
    'ito.s': 'Ito.Ito',
    'inaca.s': 'INaCa.INaCa',
    'ik1.s': 'IK1.IK1',
    'inak.s': 'INaK.INaK',
}

model_current = {
    'tnnp-2004.mmt': tnnp_2004_current,
    'fink-2008.mmt': fink_2008_current,
    'fink-2008-lei-2019-herg.mmt': fink_2008_current,
    'fink-2008-lei-2019-herg-tunable.mmt': fink_2008_current,
    'grandi-2010.mmt': grandi_2010_current,
    'ohara-2011.mmt': ohara_2011_current,
    'cipa-2017.mmt': cipa_2017_current,
    'tomek-2019.mmt': tomek_2019_current,
}


# model_conductance = {'mmt_file_name': condutance_dictionary}
# conductance_dictionary = { conductance_scaling: condutance_in_mmt }
# different model might call conductance differently...
# this is mainly for current clamp simulation (AP simulation)

tnnp_2004_conductance = {
    'ina.s': 'fast_sodium_current.g_Na',
    'ical.s': 'L_type_Ca_current.g_CaL',
    'ikr.s': 'rapid_time_dependent_potassium_current.g_Kr',
    'iks.s': 'slow_time_dependent_potassium_current.g_Ks',
    'ito.s': 'transient_outward_current.g_to',
    'inaca.s': 'sodium_calcium_exchanger_current.K_NaCa',
    'ik1.s': 'inward_rectifier_potassium_current.g_K1',
    'inak.s': 'sodium_potassium_pump_current.P_NaK',
}

fink_2008_conductance = {
    'ina.s': 'ina.g_Na',
    'ical.s': 'ical.g_CaL',
    'ikr.s': 'ikr.g_Kr_0',
    'iks.s': 'iks.g_Ks',
    'ito.s': 'ito.g_to',
    'inaca.s': 'inaca.K_NaCa',
    'ik1.s': 'ik1.g_K1_0',
    'inak.s': 'inak.P_NaK',
}

grandi_2010_conductance = {
    'ina.s': 'parameters.GNa',
    'ical.s': 'parameters.I_Ca_scale',
    'ikr.s': 'I_Kr.gkr',
    'iks.s': 'I_Ks.i_ks_scale',
    'ito.s': 'I_to.I_to_scale',
    'inaca.s': 'parameters.IbarNCX',
    'ik1.s': 'I_Ki.i_ki_scale',
    'inak.s': 'parameters.IbarNaK',
}

ohara_2011_conductance = {
    'ina.s': 'ina.GNa',
    'ical.s': 'ical.PCa',
    'ikr.s': 'ikr.GKr',
    'iks.s': 'iks.GKs',
    'ito.s': 'ito.Gto',
    'inaca.s': 'inaca.Gncx',
    'ik1.s': 'ik1.GK1',
    'inak.s': 'inak.Pnak',
}

cipa_2017_conductance = {
    'ina.s': 'INa.GNa',
    'ical.s': 'ICaL.PCa',
    'ikr.s': 'IKr.GKr',
    'iks.s': 'IKs.GKs',
    'ito.s': 'Ito.Gto',
    'inaca.s': 'INaCa_i.Gncx',
    'ik1.s': 'IK1.GK1',
    'inak.s': 'INaK.Pnak',
}

tomek_2019_conductance = {
    'ina.s': 'INa.GNa',
    'ical.s': 'ICaL.PCa',
    'ikr.s': 'IKr.GKr',
    'iks.s': 'IKs.GKs',
    'ito.s': 'Ito.Gto',
    'inaca.s': 'INaCa.Gncx',
    'ik1.s': 'IK1.GK1',
    'inak.s': 'INaK.Pnak',
}

model_conductance = {
    'tnnp-2004.mmt': tnnp_2004_conductance,
    'fink-2008.mmt': fink_2008_conductance,
    'fink-2008-lei-2019-herg.mmt': fink_2008_conductance,
    'fink-2008-lei-2019-herg-tunable.mmt': fink_2008_conductance,
    'grandi-2010.mmt': grandi_2010_conductance,
    'ohara-2011.mmt': ohara_2011_conductance,
    'cipa-2017.mmt': cipa_2017_conductance,
    'tomek-2019.mmt': tomek_2019_conductance,
}


# mainly for voltage clamp -- to clamp the ion concentration too.
default_nai = 10
default_nao = 150
default_ki = 110
default_ko = 4
default_cai = 1e-5
default_cao = 1.2
model_ion = {
    'tnnp-2004.mmt': [('sodium_dynamics.Na_i',default_nai),
                    ('sodium_dynamics.Na_o',default_nao),
                    ('potassium_dynamics.K_i',default_ki),
                    ('potassium_dynamics.K_o',default_ko),
                    ('calcium_dynamics.Ca_i',default_cai),
                    ('calcium_dynamics.Ca_o',default_cao),
                    ],
    'fink-2008.mmt': [('sodium.Na_i',default_nai),
                    ('ion.Na_o',default_nao),
                    ('potassium.K_i',default_ki),
                    ('ion.K_o',default_ko),
                    ('calcium.Ca_i',default_cai),
                    ('ion.Ca_o',default_cao),
                    ],
    'fink-2008-lei-2019-herg.mmt': [('sodium.Na_i',default_nai),
                    ('ion.Na_o',default_nao),
                    ('potassium.K_i',default_ki),
                    ('ion.K_o',default_ko),
                    ('calcium.Ca_i',default_cai),
                    ('ion.Ca_o',default_cao),
                    ],
    'fink-2008-lei-2019-herg-tunable.mmt': [('sodium.Na_i',default_nai),
                    ('ion.Na_o',default_nao),
                    ('potassium.K_i',default_ki),
                    ('ion.K_o',default_ko),
                    ('calcium.Ca_i',default_cai),
                    ('ion.Ca_o',default_cao),
                    ],
    'grandi-2010.mmt': [('Na_Concentrations.Na_i',default_nai),
                    ('parameters.Nao',default_nao),
                    ('K_Concentration.K_i',default_ki),
                    ('parameters.Ko',default_ko),
                    ('Ca_Concentrations.Ca_i',default_cai),
                    ('parameters.Cao',default_cao),
                    ],
    'ohara-2011.mmt': [('sodium.Nai',default_nai),
                    ('extra.Nao',default_nao),
                    ('potassium.Ki',default_ki),
                    ('extra.Ko',default_ko),
                    ('calcium.Cai',default_cai),
                    ('extra.Cao',default_cao),
                    ],
    'cipa-2017.mmt': [('intracellular_ions.nai',default_nai),
                    ('extracellular.nao',default_nao),
                    ('intracellular_ions.ki',default_ki),
                    ('extracellular.ko',default_ko),
                    ('intracellular_ions.cai',default_cai),
                    ('extracellular.cao',default_cao),
                    ],
    'tomek-2019.mmt': [('intracellular_ions.nai',default_nai),
                    ('extracellular.nao',default_nao),
                    ('intracellular_ions.ki',default_ki),
                    ('extracellular.ko',default_ko),
                    ('intracellular_ions.cai',default_cai),
                    ('extracellular.cao',default_cao),
                    ],
    }


# Stimulus setting
default_stim_amp = -80  # in [A/F]
model_stim_amp = {
    # 'mmt_file_name': (stim_amp, value)
    'tnnp-2004.mmt': ('membrane.stim_amplitude', default_stim_amp),
    'fink-2008.mmt': ('stimulus.stim_amplitude', default_stim_amp),
    'fink-2008-lei-2019-herg.mmt': ('stimulus.stim_amplitude', default_stim_amp),
    'fink-2008-lei-2019-herg-tunable.mmt': ('stimulus.stim_amplitude', default_stim_amp),
    'grandi-2010.mmt': ('interface.stim_amplitude', default_stim_amp),
    'ohara-2011.mmt': ('stimulus.amplitude', default_stim_amp),
    'cipa-2017.mmt': ('membrane.i_Stim_Amplitude', default_stim_amp),
    'tomek-2019.mmt': ('membrane.i_Stim_Amplitude', default_stim_amp),
    }

default_stim_setup = (1, 50, 1000, 1)
model_stim_setup = {
    # 'mmt_file_name': (stim_dur, stim_off, cycle_length, norm_stim_amp)
    'tnnp-2004.mmt': default_stim_setup,
    'fink-2008.mmt': default_stim_setup,
    'fink-2008-lei-2019-herg.mmt': default_stim_setup,
    'fink-2008-lei-2019-herg-tunable.mmt': default_stim_setup,
    'grandi-2010.mmt': default_stim_setup,
    'ohara-2011.mmt': default_stim_setup,
    'cipa-2017.mmt': default_stim_setup,
    'tomek-2019.mmt': default_stim_setup,
    }
