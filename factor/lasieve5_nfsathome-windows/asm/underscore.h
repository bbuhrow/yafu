// SMJS This file should be included at the top of asm files which use C external variables or functions.
//      It compensates for any name mangling underscores in the C variable names
#ifndef __UNDERSCORE_H__
#define __UNDERSCORE_H__

#if defined(__APPLE__)

#define mpqs3_FBk _mpqs3_FBk
#define mpqs3_FBk_inv _mpqs3_FBk_inv
#define mpqs3_FB_inv _mpqs3_FB_inv
#define mpqs3_FB_A_inv _mpqs3_FB_A_inv
#define mpqs3_FB_inv_info _mpqs3_FB_inv_info

#define mpqs_FBk _mpqs_FBk
#define mpqs_FBk_inv _mpqs_FBk_inv
#define mpqs_FB_inv _mpqs_FB_inv
#define mpqs_FB_A_inv _mpqs_FB_A_inv
#define abort _abort
#define asm3_td _asm3_td
#define asm_modinv32b _asm_modinv32b
#define asm_td _asm_td
#define asm_getbc _asm_getbc
#define get_recurrence_info _get_recurrence_info
#define modulo32 _modulo32
#define modulo64 _modulo64
#define montgomery_inv_n _montgomery_inv_n
#define montgomery_modulo_n _montgomery_modulo_n
#define mpqs_256_inv_table _mpqs_256_inv_table
#define mpqs3_FB _mpqs3_FB
#define mpqs3_FB0 _mpqs3_FB0
#define mpqs3_FB_disp _mpqs3_FB_disp
#define mpqs3_FB_log _mpqs3_FB_log
#define mpqs3_FB_mm_inv _mpqs3_FB_mm_inv
#define mpqs3_FB_np_p _mpqs3_FB_np_p
#define mpqs3_FB_np_px _mpqs3_FB_np_px
#define mpqs3_FB_start _mpqs3_FB_start
#define mpqs3_sievearray _mpqs3_sievearray
#define mpqs3_sievebegin _mpqs3_sievebegin
#define mpqs3_sievelen _mpqs3_sievelen
#define mpqs_FB _mpqs_FB
#define mpqs_FB_disp _mpqs_FB_disp
#define mpqs_FB_log _mpqs_FB_log
#define mpqs_FB_mm_inv _mpqs_FB_mm_inv
#define mpqs_FB_np_p _mpqs_FB_np_p
#define mpqs_FB_np_px _mpqs_FB_np_px
#define mpqs_FB_start _mpqs_FB_start
#define mpqs_FB_inv_info _mpqs_FB_inv_info
#define mpqs_gauss_c _mpqs_gauss_c
#define mpqs_gauss_col _mpqs_gauss_col
#define mpqs_gauss_d _mpqs_gauss_d
#define mpqs_gauss_j _mpqs_gauss_j
#define mpqs_gauss_k _mpqs_gauss_k
#define mpqs_gauss_m _mpqs_gauss_m
#define mpqs_gauss_mat _mpqs_gauss_mat
#define mpqs_gauss_n32 _mpqs_gauss_n32
#define mpqs_gauss_row _mpqs_gauss_row
#define mpqs_sievearray _mpqs_sievearray
#define mpqs_sievebegin _mpqs_sievebegin
#define mpqs_sievelen _mpqs_sievelen
#define mpqs_td_begin _mpqs_td_begin
#define mpqs_nFBk_1 _mpqs_nFBk_1
#define mpqs_nFBk _mpqs_nFBk
#define mpqs3_Adiv_all _mpqs3_Adiv_all
#define mpqs_nAdiv_total _mpqs_nAdiv_total
#define mpqs3_256_inv_table _mpqs3_256_inv_table
#define mpqs_Adiv_all _mpqs_Adiv_all
#define mpqs_nFB _mpqs_nFB
#define mpqs3_td_begin _mpqs3_td_begin
#define mpqs3_nFBk_1 _mpqs3_nFBk_1
#define mpqs3_nFBk _mpqs3_nFBk
#define mpqs3_nFB _mpqs3_nFB
#define mpqs3_nAdiv_total _mpqs3_nAdiv_total

#endif

#endif
