! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module dart_mod

! dart module to make all other module public routines/variables
! be public from this single module



!---------------------------------------------

!../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
use adaptive_inflate_mod, only : update_inflation,         do_obs_inflate,     &
          do_varying_ss_inflate,    do_single_ss_inflate,   inflate_ens,        &
          adaptive_inflate_init,    adaptive_inflate_type,                      &
                                    deterministic_inflate,  solve_quadratic,    &
          log_inflation_info,       get_minmax_task_zero,   mean_from_restart,  &
          sd_from_restart,                                                      &
          output_inf_restart,       get_inflate_mean,       get_inflate_sd,     &
          get_is_prior,             get_is_posterior,       do_ss_inflate,      &
          set_inflation_mean_copy,  set_inflation_sd_copy,  get_inflation_mean_copy, &
          get_inflation_sd_copy,    do_rtps_inflate,        validate_inflate_options, &
          print_inflation_restart_filename, &
          PRIOR_INF, POSTERIOR_INF, NO_INFLATION, OBS_INFLATION, VARYING_SS_INFLATION, &
          SINGLE_SS_INFLATION, RELAXATION_TO_PRIOR_SPREAD, ENHANCED_SS_INFLATION

!../assimilation_code/modules/assimilation/assim_model_mod.f90
use assim_model_mod, only : static_init_assim_model, &
          end_assim_model, &
          get_initial_condition, &
          get_closest_state_time_to, &
          get_state_meta_data, &
          get_model_time_step, &
          get_model_size,  &
          adv_1step, &
          interpolate, &
          pert_model_copies, &
          get_close_obs, &
          get_close_state, &
          convert_vertical_obs, &
          convert_vertical_state, &
          read_model_time, &
          write_model_time

!../assimilation_code/modules/assimilation/assim_tools_mod.f90
use assim_tools_mod, only : filter_assim, &
          set_assim_tools_trace, &
          test_state_copies, &
          update_ens_from_weights  ! Jeff thinks this routine is in the wild.

!../assimilation_code/modules/assimilation/cov_cutoff_mod.f90
use cov_cutoff_mod, only : comp_cov_factor

!../assimilation_code/modules/io/dart_time_io_mod.f90
use dart_time_io_mod, only : read_model_time, write_model_time

!../assimilation_code/modules/io/direct_netcdf_mod.f90
use direct_netcdf_mod, only : read_transpose,            &
          transpose_write,           &
          initialize_single_file_io, &
          finalize_single_file_io,   &
          read_single_file,          &
          write_single_file,         &
          write_augmented_state,     &
          read_variables,            &
          nc_get_num_times


!../assimilation_code/modules/utilities/distributed_state_mod.f90
use distributed_state_mod, only : get_state_array, get_state, create_state_window, &
          free_state_window, create_mean_window, free_mean_window

!../assimilation_code/modules/utilities/ensemble_manager_mod.f90
use ensemble_manager_mod, only : copies_in_window, mean_row, set_num_extra_copies,   &
          get_allow_transpose
use ensemble_manager_mod, only : init_ensemble_manager,      end_ensemble_manager,     get_ensemble_time,          &
          ensemble_type,              duplicate_ens,            get_var_owner_index,        &
          get_my_num_copies,          get_my_copies,            get_my_num_vars,            &
          get_my_vars,                compute_copy_mean,        compute_copy_mean_sd,       &
          get_copy,                   put_copy,                 all_vars_to_all_copies,     &
          all_copies_to_all_vars,     allocate_vars,            deallocate_vars,            &
          compute_copy_mean_var,      get_copy_owner_index,     set_ensemble_time,          &
          broadcast_copy,             prepare_to_write_to_vars, prepare_to_write_to_copies, &
          prepare_to_read_from_vars,  prepare_to_read_from_copies, prepare_to_update_vars,  &
          prepare_to_update_copies,   print_ens_handle,         set_current_time,           &
          map_task_to_pe,             map_pe_to_task,           get_current_time,           &
          allocate_single_copy,       put_single_copy,          get_single_copy,            &
          deallocate_single_copy

!../assimilation_code/modules/assimilation/filter_mod.f90
use filter_mod, only : filter_sync_keys_time, &
          filter_set_initial_time, &
          filter_main

!../assimilation_code/modules/observations/forward_operator_mod.f90
use forward_operator_mod, only : get_obs_ens_distrib_state, get_expected_obs_distrib_state

!../assimilation_code/modules/io/io_filenames_mod.f90
use io_filenames_mod, only : io_filenames_init, &
          io_filenames_finalize, &
          file_info_type, &
          netcdf_file_type, &
          stage_metadata_type, &
          set_file_metadata, &
          set_member_file_metadata, &
          set_io_copy_flag, &
          assert_file_info_initialized, &
          assert_restart_names_initialized, &
          file_info_dump, &
          combine_file_info, &
          check_file_info_variable_shape
use io_filenames_mod, only : get_restart_filename, &
          get_single_file, &
          get_cycling, &
          get_file_description, &
          get_copy_name, &
          get_stage_metadata, &
          single_file_initialized, &
          inherit_copy_units, &
          copy_is_clamped, &
          force_copy_back, &
          noutput_state_variables
use io_filenames_mod, only : query_read_copy, &
          query_write_copy, &
          query_copy_present
use io_filenames_mod, only : READ_COPY, &
          WRITE_COPY, &
          READ_WRITE_COPY, &
          NO_IO, &
          COPY_NOT_PRESENT

!../assimilation_code/modules/utilities/netcdf_utilities_mod.f90
use netcdf_utilities_mod, only : nc_check,                       &
          nc_add_global_attribute,        &
          nc_get_global_attribute,        &
          nc_add_attribute_to_variable,   &
          nc_get_attribute_from_variable, &
          nc_define_dimension,            &
          nc_define_unlimited_dimension,  &
          nc_get_dimension_size,          &
          nc_define_character_variable,   &
          nc_define_integer_variable,     &
          nc_define_real_variable,        &
          nc_define_double_variable,      &
          nc_define_integer_scalar,       &
          nc_define_real_scalar,          &
          nc_define_double_scalar,        &
          nc_global_attribute_exists,     &
          nc_variable_attribute_exists,   &
          nc_dimension_exists,            &
          nc_variable_exists,             &
          nc_put_variable,                &
          nc_get_variable,                &
          nc_get_variable_info,           &
          nc_add_global_creation_time,    &
          nc_get_variable_num_dimensions, &
          nc_get_variable_dimension_names, &
          nc_get_variable_size,           &
          nc_open_file_readonly,          &
          nc_open_file_readwrite,         &
          nc_create_file,                 &
          nc_close_file,                  &
          nc_begin_define_mode,           &
          nc_end_define_mode,             &
          nc_synchronize_file,            &
          NF90_MAX_NAME, NF90_MAX_VAR_DIMS

!../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, get_dart_mpi_comm, datasize, send_minmax_to,     &
          get_from_fwd, get_from_mean, broadcast_minmax, broadcast_flag,     &
          start_mpi_timer, read_mpi_timer, send_sum_to, get_global_max,      &
          all_reduce_min_max  ! deprecated, replace by broadcast_minmax

!../observations/forward_operators/obs_def_utilities_mod.f90
use obs_def_utilities_mod, only : track_status, set_debug_fwd_op

!../assimilation_code/modules/utilities/obs_impact_mod.f90
use obs_impact_mod, only : create_impact_table,   &
          allocate_impact_table, &
          read_impact_table,     &
          free_impact_table

!../assimilation_code/modules/assimilation/obs_model_mod.f90
use obs_model_mod, only : move_ahead, advance_state, set_obs_model_trace, have_members

!../assimilation_code/modules/observations/obs_sequence_mod.f90
use obs_sequence_mod, only : obs_sequence_type, init_obs_sequence, interactive_obs_sequence, &
   get_num_copies, get_num_qc, get_num_obs, get_max_num_obs, &
   get_copy_meta_data, get_qc_meta_data, get_next_obs, get_prev_obs, &
   insert_obs_in_seq, delete_obs_from_seq, set_copy_meta_data, &
   set_qc_meta_data, get_first_obs, get_last_obs, add_copies, add_qc, &
   write_obs_seq, read_obs_seq, set_obs, append_obs_to_seq, &
   get_obs_from_key, get_obs_time_range, get_time_range_keys, &
   get_num_times, get_num_key_range, operator(==), operator(/=), &
   static_init_obs_sequence, destroy_obs_sequence, read_obs_seq_header, &
   delete_seq_head, delete_seq_tail, &
   get_next_obs_from_key, get_prev_obs_from_key, delete_obs_by_typelist, &
   select_obs_by_location, delete_obs_by_qc, delete_obs_by_copy, &
   print_obs_seq_summary, validate_obs_seq_time
use obs_sequence_mod, only : obs_type, init_obs, destroy_obs, get_obs_def, set_obs_def, &
   get_obs_values, set_obs_values, replace_obs_values, get_qc, set_qc, &  
   read_obs, write_obs, replace_qc, interactive_obs, copy_obs, assignment(=), &
   get_obs_key, copy_partial_obs, print_obs
use obs_sequence_mod, only : obs_cov_type

!../assimilation_code/modules/utilities/options_mod.f90
use options_mod, only : get_missing_ok_status, set_missing_ok_status

!../assimilation_code/modules/utilities/parse_args_mod.f90
use parse_args_mod, only : get_args_from_string, &
          get_name_val_pairs_from_string, &
          get_next_arg

!../assimilation_code/modules/assimilation/quality_control_mod.f90
use quality_control_mod, only : initialize_qc, input_qc_ok, get_dart_qc, check_outlier_threshold, &
          good_dart_qc, set_input_qc, &
          DARTQC_ASSIM_GOOD_FOP, DARTQC_EVAL_GOOD_FOP, &
          DARTQC_ASSIM_FAILED_POST_FOP, DARTQC_EVAL_FAILED_POST_FOP, &
          DARTQC_FAILED_FOP, DARTQC_NOT_IN_NAMELIST, &
          DARTQC_BAD_INCOMING_QC, DARTQC_FAILED_OUTLIER_TEST, &
          DARTQC_FAILED_VERT_CONVERT

!../assimilation_code/modules/utilities/random_seq_mod.f90
use random_seq_mod, only : random_seq_type, &
          init_random_seq, &
          random_uniform, &
          random_gaussian, &
          several_random_gaussians, &
          twod_gaussians, &
          random_gamma, & 
          random_inverse_gamma, & 
          random_exponential

!../assimilation_code/modules/assimilation/reg_factor_mod.f90
use reg_factor_mod, only : comp_reg_factor

!../assimilation_code/modules/assimilation/sampling_error_correction_mod.f90
use sampling_error_correction_mod, only : get_sampling_error_table_size, &
          read_sampling_error_correction

!../assimilation_code/modules/assimilation/smoother_mod.f90
use smoother_mod, only : smoother_read_restart, advance_smoother,                     &
   smoother_gen_copy_meta_data, smoother_write_restart, init_smoother, &
   do_smoothing, smoother_mean_spread, smoother_assim,                 &
   smoother_ss_diagnostics,            &
   smoother_end, set_smoother_trace

!../assimilation_code/modules/utilities/sort_mod.f90
use sort_mod, only : sort, index_sort, insertion_sort, index_insertion_sort
!use sort_mod, only : simple_sort, simple_index_sort  

!../assimilation_code/modules/io/state_structure_mod.f90
use state_structure_mod, only : add_domain,                 &
          get_domain_size,            &
          get_num_domains,            &
          get_variable_size,          &
          get_variable_name,          &
          get_kind_string,            &
          get_kind_index,             &
          get_varid_from_varname,     & 
          get_varid_from_kind,        &
          get_varids_from_kind,       &
          get_num_variables,          &
          get_num_dims,               &
          get_dim_lengths,            &
          get_dim_length,             &
          get_dim_name,               &
          get_io_num_dims,            &
          get_io_dim_ids,             &
          get_io_dim_lengths,         &
          get_io_num_unique_dims,     &
          get_io_unique_dim_name,     &
          get_io_unique_dim_length,   &
          add_time_unlimited,         &
          get_unlimited_dimid,        &
          set_var_id,                 &
          get_io_clamping_maxval,     &
          get_io_clamping_minval,     &
          do_io_clamping,             &
          do_io_update,               &
          get_index_start,            &
          get_index_end,              &
          get_sum_variables,          &
          get_sum_variables_below,    &
          get_model_variable_indices, &
          get_dart_vector_index,      &
          get_num_varids_from_kind,   &
          get_xtype,                  &
          get_units,                  &
          get_long_name,              &
          get_short_name,             &
          get_has_missing_value,      &
          get_FillValue,              &
          get_missing_value,          &
          get_add_offset,             &
          get_scale_factor,           &
          set_dart_kinds,             &
          set_clamping,               &
          set_update_list,            &
          add_dimension_to_variable,  &
          finished_adding_domain,     &
          state_structure_info
use state_structure_mod, only : create_diagnostic_structure, &
          end_diagnostic_structure

!../assimilation_code/modules/io/state_vector_io_mod.f90
use state_vector_io_mod, only : state_vector_io_init, &
          read_state, &
          write_state
use state_vector_io_mod, only : set_stage_to_write, &
          get_stage_to_write

!../assimilation_code/modules/utilities/time_manager_mod.f90
use time_manager_mod, only : time_type

! Operators defined on time_type
use time_manager_mod, only : operator(+),  operator(-),   operator(*),   operator(/),  &
          operator(>),  operator(>=),  operator(==),  operator(/=), &
          operator(<),  operator(<=),  operator(//)
use time_manager_mod, only : set_time, set_time_missing, increment_time, decrement_time, get_time
use time_manager_mod, only : interval_alarm, repeat_alarm, generate_seed, time_index_sort
use time_manager_mod, only : THIRTY_DAY_MONTHS,    JULIAN,    GREGORIAN,  NOLEAP,   NO_CALENDAR, &
          GREGORIAN_MARS,   SOLAR_MARS
use time_manager_mod, only : set_calendar_type, get_calendar_type, get_calendar_string
use time_manager_mod, only : set_date,       set_date_gregorian,         set_date_julian, &
                          set_date_thirty,            set_date_no_leap, &
                          set_date_gregorian_mars, set_date_solar_mars
use time_manager_mod, only : get_date,       get_date_gregorian,         get_date_julian, &
                          get_date_thirty,            get_date_no_leap, &
                          get_date_gregorian_mars, get_date_solar_mars
use time_manager_mod, only : increment_date, increment_gregorian,        increment_julian, &
                          increment_thirty,           increment_no_leap, &
                          increment_gregorian_mars, increment_solar_mars
use time_manager_mod, only : decrement_date, decrement_gregorian,        decrement_julian, &
                          decrement_thirty,           decrement_no_leap, &
                          decrement_gregorian_mars, decrement_solar_mars
use time_manager_mod, only : days_in_month,  days_in_month_gregorian,    days_in_month_julian, &
                          days_in_month_no_leap,      days_in_month_thirty, &
                          days_in_month_gregorian_mars, days_in_month_solar_mars
use time_manager_mod, only : leap_year,      leap_year_gregorian,        leap_year_julian, &
                          leap_year_no_leap,          leap_year_thirty, &
                          leap_year_gregorian_mars, leap_year_solar_mars
use time_manager_mod, only : length_of_year, length_of_year_thirty,      length_of_year_julian, &
                          length_of_year_gregorian,   length_of_year_no_leap, &
                          length_of_year_gregorian_mars, length_of_year_solar_mars
use time_manager_mod, only : days_in_year,   days_in_year_thirty,        days_in_year_julian, &
                          days_in_year_gregorian,     days_in_year_no_leap, &
                          days_in_year_gregorian_mars, days_in_year_solar_mars
use time_manager_mod, only : month_name

use time_manager_mod, only : julian_day
use time_manager_mod, only : time_manager_init, print_time, print_date
use time_manager_mod, only : write_time, read_time, interactive_time

!../assimilation_code/modules/utilities/types_mod.f90
use types_mod, only : i2, i4, i8, r4, c4, r8, c8, digits12
use types_mod, only : PI, DEG2RAD, RAD2DEG, MISSING_R4, MISSING_R8
use types_mod, only : MISSING_I, MISSING_I8, MISSING_DATA
use types_mod, only : SECPERDAY
use types_mod, only : t_kelvin, es_alpha, es_beta, es_gamma, gas_constant_v, gas_constant
use types_mod, only : L_over_Rv, ps0, earth_radius, gravity
use types_mod, only : metadatalength, obstypelength, varnamelength, vtablenamelength
use types_mod, only : MAX_NUM_DOMS, MAX_FILES

!---------------------------------------------

implicit none
private

! Operators.  these need to be made public only once regardless of
! the underlying type.  only list each type once, and only list ones
! where at least one module supplies it.
public :: operator(+),  operator(-),   operator(*),   operator(/),  &
          operator(>),  operator(>=),  operator(==),  operator(/=), &
          operator(<),  operator(<=),  operator(//)

! ../assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
public :: update_inflation,                                 do_obs_inflate,     &
          do_varying_ss_inflate,    do_single_ss_inflate,   inflate_ens,        &
          adaptive_inflate_init,    adaptive_inflate_type,                      &
                                    deterministic_inflate,  solve_quadratic,    &
          log_inflation_info,       get_minmax_task_zero,   mean_from_restart,  &
          sd_from_restart,                                                      &
          output_inf_restart,       get_inflate_mean,       get_inflate_sd,     &
          get_is_prior,             get_is_posterior,       do_ss_inflate,      &
          set_inflation_mean_copy,  set_inflation_sd_copy,  get_inflation_mean_copy, &
          get_inflation_sd_copy,    do_rtps_inflate,        validate_inflate_options, &
          print_inflation_restart_filename, &
          PRIOR_INF, POSTERIOR_INF, NO_INFLATION, OBS_INFLATION, VARYING_SS_INFLATION, &
          SINGLE_SS_INFLATION, RELAXATION_TO_PRIOR_SPREAD, ENHANCED_SS_INFLATION

!../assimilation_code/modules/assimilation/assim_model_mod.f90
public :: static_init_assim_model, &
          end_assim_model

!../assimilation_code/modules/assimilation/assim_tools_mod.f90
public :: filter_assim, &
          set_assim_tools_trace, &
          test_state_copies, &
          update_ens_from_weights  ! Jeff thinks this routine is in the wild.

!../assimilation_code/modules/assimilation/cov_cutoff_mod.f90
public :: comp_cov_factor

!../assimilation_code/modules/io/dart_time_io_mod.f90
public :: read_model_time, write_model_time

!../assimilation_code/modules/io/direct_netcdf_mod.f90
public :: read_transpose,            &
          transpose_write,           &
          initialize_single_file_io, &
          finalize_single_file_io,   &
          read_single_file,          &
          write_single_file,         &
          write_augmented_state,     &
          read_variables,            &
          nc_get_num_times


!../assimilation_code/modules/utilities/distributed_state_mod.f90
public :: get_state_array, get_state, create_state_window, &
          free_state_window, create_mean_window, free_mean_window

!../assimilation_code/modules/utilities/ensemble_manager_mod.f90
public :: copies_in_window, mean_row, set_num_extra_copies,   &
          get_allow_transpose
public :: init_ensemble_manager,      end_ensemble_manager,     get_ensemble_time,         &
          ensemble_type,              duplicate_ens,            get_var_owner_index,        &
          get_my_num_copies,          get_my_copies,            get_my_num_vars,            &
          get_my_vars,                compute_copy_mean,        compute_copy_mean_sd,       &
          get_copy,                   put_copy,                 all_vars_to_all_copies,     &
          all_copies_to_all_vars,     allocate_vars,            deallocate_vars,            &
          compute_copy_mean_var,      get_copy_owner_index,     set_ensemble_time,          &
          broadcast_copy,             prepare_to_write_to_vars, prepare_to_write_to_copies, &
          prepare_to_read_from_vars,  prepare_to_read_from_copies, prepare_to_update_vars,  &
          prepare_to_update_copies,   print_ens_handle,         set_current_time,           &
          map_task_to_pe,             map_pe_to_task,           get_current_time,           &
          allocate_single_copy,       put_single_copy,          get_single_copy,            &
          deallocate_single_copy

!../assimilation_code/modules/assimilation/filter_mod.f90
public :: filter_sync_keys_time, &
          filter_set_initial_time, &
          filter_main

!../assimilation_code/modules/observations/forward_operator_mod.f90
public :: get_obs_ens_distrib_state, get_expected_obs_distrib_state

!../assimilation_code/modules/io/io_filenames_mod.f90
public :: io_filenames_init, &
          io_filenames_finalize, &
          file_info_type, &
          netcdf_file_type, &
          stage_metadata_type, &
          set_file_metadata, &
          set_member_file_metadata, &
          set_io_copy_flag, &
          assert_file_info_initialized, &
          assert_restart_names_initialized, &
          file_info_dump, &
          combine_file_info, &
          check_file_info_variable_shape
public :: get_restart_filename, &
          get_single_file, &
          get_cycling, &
          get_file_description, &
          get_copy_name, &
          get_stage_metadata, &
          single_file_initialized, &
          inherit_copy_units, &
          copy_is_clamped, &
          force_copy_back, &
          noutput_state_variables
public :: query_read_copy, &
          query_write_copy, &
          query_copy_present
public :: READ_COPY, &
          WRITE_COPY, &
          READ_WRITE_COPY, &
          NO_IO, &
          COPY_NOT_PRESENT

!../assimilation_code/modules/utilities/netcdf_utilities_mod.f90
public :: nc_check,                       &
          nc_add_global_attribute,        &
          nc_get_global_attribute,        &
          nc_add_attribute_to_variable,   &
          nc_get_attribute_from_variable, &
          nc_define_dimension,            &
          nc_define_unlimited_dimension,  &
          nc_get_dimension_size,          &
          nc_define_character_variable,   &
          nc_define_integer_variable,     &
          nc_define_real_variable,        &
          nc_define_double_variable,      &
          nc_define_integer_scalar,       &
          nc_define_real_scalar,          &
          nc_define_double_scalar,        &
          nc_global_attribute_exists,     &
          nc_variable_attribute_exists,   &
          nc_dimension_exists,            &
          nc_variable_exists,             &
          nc_put_variable,                &
          nc_get_variable,                &
          nc_get_variable_info,           &
          nc_add_global_creation_time,    &
          nc_get_variable_num_dimensions, &
          nc_get_variable_dimension_names, &
          nc_get_variable_size,           &
          nc_open_file_readonly,          &
          nc_open_file_readwrite,         &
          nc_create_file,                 &
          nc_close_file,                  &
          nc_begin_define_mode,           &
          nc_end_define_mode,             &
          nc_synchronize_file,            &
          NF90_MAX_NAME, NF90_MAX_VAR_DIMS

!../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
public :: initialize_mpi_utilities, finalize_mpi_utilities,                  &
          task_count, my_task_id, block_task, restart_task,                  &
          task_sync, array_broadcast, send_to, receive_from, iam_task0,      &
          broadcast_send, broadcast_recv, shell_execute, sleep_seconds,      &
          sum_across_tasks, get_dart_mpi_comm, datasize, send_minmax_to,     &
          get_from_fwd, get_from_mean, broadcast_minmax, broadcast_flag,     &
          start_mpi_timer, read_mpi_timer, send_sum_to, get_global_max,      &
          all_reduce_min_max  ! deprecated, replace by broadcast_minmax

!../observations/forward_operators/obs_def_utilities_mod.f90
public :: track_status, set_debug_fwd_op

!../assimilation_code/modules/utilities/obs_impact_mod.f90
public :: create_impact_table,   &
          allocate_impact_table, &
          read_impact_table,     &
          free_impact_table

!../assimilation_code/modules/assimilation/obs_model_mod.f90
public :: move_ahead, advance_state, set_obs_model_trace, have_members

!../assimilation_code/modules/observations/obs_sequence_mod.f90
public :: obs_sequence_type, init_obs_sequence, interactive_obs_sequence, &
   get_num_copies, get_num_qc, get_num_obs, get_max_num_obs, &
   get_copy_meta_data, get_qc_meta_data, get_next_obs, get_prev_obs, &
   insert_obs_in_seq, delete_obs_from_seq, set_copy_meta_data, &
   set_qc_meta_data, get_first_obs, get_last_obs, add_copies, add_qc, &
   write_obs_seq, read_obs_seq, set_obs, append_obs_to_seq, &
   get_obs_from_key, get_obs_time_range, get_time_range_keys, &
   get_num_times, get_num_key_range, & !!! operator(==), operator(/=), &
   static_init_obs_sequence, destroy_obs_sequence, read_obs_seq_header, &
   delete_seq_head, delete_seq_tail, &
   get_next_obs_from_key, get_prev_obs_from_key, delete_obs_by_typelist, &
   select_obs_by_location, delete_obs_by_qc, delete_obs_by_copy, &
   print_obs_seq_summary, validate_obs_seq_time
public :: obs_type, init_obs, destroy_obs, get_obs_def, set_obs_def, &
   get_obs_values, set_obs_values, replace_obs_values, get_qc, set_qc, &  
   read_obs, write_obs, replace_qc, interactive_obs, copy_obs, assignment(=), &
   get_obs_key, copy_partial_obs, print_obs
public :: obs_cov_type

!../assimilation_code/modules/utilities/options_mod.f90
public :: get_missing_ok_status, set_missing_ok_status

!../assimilation_code/modules/utilities/parse_args_mod.f90
public :: get_args_from_string, &
          get_name_val_pairs_from_string, &
          get_next_arg

!../assimilation_code/modules/assimilation/quality_control_mod.f90
public :: initialize_qc, input_qc_ok, get_dart_qc, check_outlier_threshold, &
          good_dart_qc, set_input_qc, &
          DARTQC_ASSIM_GOOD_FOP, DARTQC_EVAL_GOOD_FOP, &
          DARTQC_ASSIM_FAILED_POST_FOP, DARTQC_EVAL_FAILED_POST_FOP, &
          DARTQC_FAILED_FOP, DARTQC_NOT_IN_NAMELIST, &
          DARTQC_BAD_INCOMING_QC, DARTQC_FAILED_OUTLIER_TEST, &
          DARTQC_FAILED_VERT_CONVERT

!../assimilation_code/modules/utilities/random_seq_mod.f90
public :: random_seq_type, &
          init_random_seq, &
          random_uniform, &
          random_gaussian, &
          several_random_gaussians, &
          twod_gaussians, &
          random_gamma, & 
          random_inverse_gamma, & 
          random_exponential

!../assimilation_code/modules/assimilation/reg_factor_mod.f90
public :: comp_reg_factor

!../assimilation_code/modules/assimilation/sampling_error_correction_mod.f90
public :: get_sampling_error_table_size, &
          read_sampling_error_correction

!../assimilation_code/modules/assimilation/smoother_mod.f90
public :: smoother_read_restart, advance_smoother,                     &
   smoother_gen_copy_meta_data, smoother_write_restart, init_smoother, &
   do_smoothing, smoother_mean_spread, smoother_assim,                 &
   smoother_ss_diagnostics,            &
   smoother_end, set_smoother_trace

!../assimilation_code/modules/utilities/sort_mod.f90
public :: sort, index_sort, insertion_sort, index_insertion_sort

!../assimilation_code/modules/io/state_structure_mod.f90
public :: add_domain,                 &
          get_domain_size,            &
          get_num_domains,            &
          get_variable_size,          &
          get_variable_name,          &
          get_kind_string,            &
          get_kind_index,             &
          get_varid_from_varname,     & 
          get_varid_from_kind,        &
          get_varids_from_kind,       &
          get_num_variables,          &
          get_num_dims,               &
          get_dim_lengths,            &
          get_dim_length,             &
          get_dim_name,               &
          get_io_num_dims,            &
          get_io_dim_ids,             &
          get_io_dim_lengths,         &
          get_io_num_unique_dims,     &
          get_io_unique_dim_name,     &
          get_io_unique_dim_length,   &
          add_time_unlimited,         &
          get_unlimited_dimid,        &
          set_var_id,                 &
          get_io_clamping_maxval,     &
          get_io_clamping_minval,     &
          do_io_clamping,             &
          do_io_update,               &
          get_index_start,            &
          get_index_end,              &
          get_sum_variables,          &
          get_sum_variables_below,    &
          get_model_variable_indices, &
          get_dart_vector_index,      &
          get_num_varids_from_kind,   &
          get_xtype,                  &
          get_units,                  &
          get_long_name,              &
          get_short_name,             &
          get_has_missing_value,      &
          get_FillValue,              &
          get_missing_value,          &
          get_add_offset,             &
          get_scale_factor,           &
          set_dart_kinds,             &
          set_clamping,               &
          set_update_list,            &
          add_dimension_to_variable,  &
          finished_adding_domain,     &
          state_structure_info
public :: create_diagnostic_structure, &
          end_diagnostic_structure

!../assimilation_code/modules/io/state_vector_io_mod.f90
public :: state_vector_io_init, &
          read_state, &
          write_state
public :: set_stage_to_write, &
          get_stage_to_write

!../assimilation_code/modules/utilities/time_manager_mod.f90
public :: time_type

public :: set_time, set_time_missing, increment_time, decrement_time, get_time
public :: interval_alarm, repeat_alarm, generate_seed, time_index_sort
public :: THIRTY_DAY_MONTHS,    JULIAN,    GREGORIAN,  NOLEAP,   NO_CALENDAR, &
          GREGORIAN_MARS,   SOLAR_MARS
public :: set_calendar_type, get_calendar_type, get_calendar_string
public :: set_date,       set_date_gregorian,         set_date_julian, &
                          set_date_thirty,            set_date_no_leap, &
                          set_date_gregorian_mars, set_date_solar_mars
public :: get_date,       get_date_gregorian,         get_date_julian, &
                          get_date_thirty,            get_date_no_leap, &
                          get_date_gregorian_mars, get_date_solar_mars
public :: increment_date, increment_gregorian,        increment_julian, &
                          increment_thirty,           increment_no_leap, &
                          increment_gregorian_mars, increment_solar_mars
public :: decrement_date, decrement_gregorian,        decrement_julian, &
                          decrement_thirty,           decrement_no_leap, &
                          decrement_gregorian_mars, decrement_solar_mars
public :: days_in_month,  days_in_month_gregorian,    days_in_month_julian, &
                          days_in_month_no_leap,      days_in_month_thirty, &
                          days_in_month_gregorian_mars, days_in_month_solar_mars
public :: leap_year,      leap_year_gregorian,        leap_year_julian, &
                          leap_year_no_leap,          leap_year_thirty, &
                          leap_year_gregorian_mars, leap_year_solar_mars
public :: length_of_year, length_of_year_thirty,      length_of_year_julian, &
                          length_of_year_gregorian,   length_of_year_no_leap, &
                          length_of_year_gregorian_mars, length_of_year_solar_mars
public :: days_in_year,   days_in_year_thirty,        days_in_year_julian, &
                          days_in_year_gregorian,     days_in_year_no_leap, &
                          days_in_year_gregorian_mars, days_in_year_solar_mars
public :: month_name

public :: julian_day
public :: time_manager_init, print_time, print_date
public :: write_time, read_time, interactive_time

!../assimilation_code/modules/utilities/types_mod.f90
public :: i2, i4, i8, r4, c4, r8, c8, digits12
public :: PI, DEG2RAD, RAD2DEG, MISSING_R4, MISSING_R8
public :: MISSING_I, MISSING_I8, MISSING_DATA
public :: SECPERDAY
public :: t_kelvin, es_alpha, es_beta, es_gamma, gas_constant_v, gas_constant
public :: L_over_Rv, ps0, earth_radius, gravity
public :: metadatalength, obstypelength, varnamelength, vtablenamelength
public :: MAX_NUM_DOMS, MAX_FILES

!---------------------------------------------

end module dart_mod
