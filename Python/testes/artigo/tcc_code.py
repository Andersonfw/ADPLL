for k in range(1 , TIME):
        K += FCW  # reference phase accumulator
        RR_k += FCW
        t_R = k * FREF_edge
        while t_CKV[n] < t_R:
            n += 1
            RV_n = n  # variable phase accumulator
            jitter.append(np.random.randn() * Jt_noise)
            wander.append(np.random.randn() * Wt_noise)
            t_CKV.append(t_CKV[n - 1] + T0 + jitter[n] + wander[n] - jitter[n - 1])  # - TDEV_I
        if trk_bank_calib:
            error_fractional[k] = TDC(t_R , t_CKV[n] ,1 / (np.sum(freq_array)))
        RV_k = RV_n  # variable phase accumulator
        phase_error[k] = (RR_k - RV_k + error_fractional[k])  # Phase detector
        ##################### PVT MODE #################################################
        if not pvt_bank_calib:
            NTW[k] = (int(phase_error[k]) * Kp_PVT)  # calcula o novo valor de NTW como inteiro
            OTW[k] = NTW[k] * (FREF / FREQ_RES_PVT)  # gain normalization of TRK mode
            OTW_pvt = OTW[k]
            if diff_freq(f_CKV) <= FREQ_RES_PVT:
                count += 1
                if count == 10:
                    pvt_bank_calib = True
                    count = 0
            else:
                count = 0
        ##################### ACQUISITION MODE #########################################
        elif not acq_bank_calib:
            NTW[k] = (int(phase_error[k]) * Kp_ACQ)  # calcula o novo valor de NTW como inteiro
            OTW[k] = NTW[k] * (FREF / FREQ_RES_ACQ) # gain normalization of TRK mode
            if OTW[k] > 2**ACQ_NB:
                OTW[k] = 2**ACQ_NB
            OTW_acq = OTW[k]  # ajusta o novo valor de controle dos capacitores do ACQ bank
            if diff_freq(f_CKV) <= FREQ_RES_ACQ:
                count += 1
                if count == 10:
                    acq_bank_calib = True
                    trk_bank_calib = True
                    OTW_acq = OTW[k - 1]
                    OFFSET_ERROR_ACQ = int(phase_error[k])
            else:
                count = 0
        ##################### TREKING MODE ################################################
        elif trk_bank_calib:
            phase_error[k] = abs(phase_error[k]) - OFFSET_ERROR_ACQ
            fractional_error_trk.append(phase_error[k])
            phase_error[k] = IRR_filter(k , phase_error[k])  # aplica o filtro IRR
            NTW[k] = (phase_error[k]) * Kp_TRK - Kp_TRK * (phase_error[k - 1]) + Ki_TRK * (phase_error[k - 1]) + NTW[k - 1]  # calcula o novo valor de NTW
            OTW[k] = NTW[k] * (FREF / FREQ_RES_TRK) # gain normalization of TRK mode
            if OTW[k] > 2**TRK_NB_I:
                OTW[k] = 2**TRK_NB_I
            elif OTW[k] < 0:
                OTW[k] = 0
            OTW_trk = OTW[k]  # calcula o novo valor de NTW como inteiro
        #######################################################################################################################
        f_CKV = SET_DCO(OTW_pvt , OTW_acq , OTW_trk )
        last_To = T0
        T0 = 1 / f_CKV
        freqs[k] = f_CKV/DIVISION_OUTPUT  # insere o valor de frequência ajustado no index k