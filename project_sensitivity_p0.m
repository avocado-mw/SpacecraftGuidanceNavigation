%-------------------------------------------------------------------------%
%    project_sensitivity_p0
%
%  DESCRIPTION:
%     Sensitivity study for the a priori station-coordinate uncertainty of
%     stations 337 and 394.  Runs batch and sequential filters for several
%     P0 scalings and compares post-fit RMS and position covariance.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
clearvars; close all; clc

S = project_setup();
S.makePlots = false;

stationVars = [ 1e4 , 1e6 , 1e8 ];
nCases      = numel( stationVars );

batch.rho_rms_post    = zeros( nCases , 1 );
batch.rhodot_rms_post = zeros( nCases , 1 );
batch.sigma_rms       = zeros( nCases , 1 );
batch.trace           = zeros( nCases , 1 );

seq.rho_rms_post    = zeros( nCases , 1 );
seq.rhodot_rms_post = zeros( nCases , 1 );
seq.sigma_rms       = zeros( nCases , 1 );
seq.trace           = zeros( nCases , 1 );

fprintf( 'P0 sensitivity: station 337/394 coordinate variance sweep\n' );
fprintf( 'Station 101 variance fixed at 1e-10 m^2\n\n' );
fprintf( '%10s | %12s %12s %12s %12s | %12s %12s %12s %12s\n' , ...
    'Var[s2,s3]' , 'Batch rho' , 'Batch rhodot' , 'Batch sig' , 'Batch tr' , ...
    'Seq rho' , 'Seq rhodot' , 'Seq sig' , 'Seq tr' );
fprintf( '%s\n' , repmat( '-' , 1 , 110 ) );

for i = 1:nCases

    P0_bar = build_P0_bar( stationVars(i) );

    batchOut = run_batch_filter( S , P0_bar );
    seqOut   = run_sequential_filter( S , P0_bar );

    batch.rho_rms_post(i)    = batchOut.rho_rms_post;
    batch.rhodot_rms_post(i) = batchOut.rhodot_rms_post;
    batch.sigma_rms(i)       = batchOut.pos_metrics.sigma_rms;
    batch.trace(i)           = batchOut.pos_metrics.trace;

    seq.rho_rms_post(i)    = seqOut.rho_rms_post;
    seq.rhodot_rms_post(i) = seqOut.rhodot_rms_post;
    seq.sigma_rms(i)       = seqOut.pos_metrics.sigma_rms;
    seq.trace(i)           = seqOut.pos_metrics.trace;

    fprintf( '%10.0e | %12.4f %12.6f %12.3f %12.3e | %12.4f %12.6f %12.3f %12.3e\n' , ...
        stationVars(i) , ...
        batchOut.rho_rms_post , batchOut.rhodot_rms_post , batchOut.pos_metrics.sigma_rms , batchOut.pos_metrics.trace , ...
        seqOut.rho_rms_post , seqOut.rhodot_rms_post , seqOut.pos_metrics.sigma_rms , seqOut.pos_metrics.trace );

end

labels = arrayfun( @(v) sprintf( 'P0=%.0e' , v ) , stationVars , 'UniformOutput' , false );

figure( 'Name' , 'P0 Sensitivity: Post-fit Range RMS' );
bar( [ batch.rho_rms_post , seq.rho_rms_post ] );
set( gca , 'XTickLabel' , labels );
ylabel( 'post-fit range RMS [m]' );
legend( 'Batch' , 'Sequential' , 'Location' , 'best' );
grid on
title( 'Post-fit Range RMS vs Station P0 Variance' );

figure( 'Name' , 'P0 Sensitivity: Position Covariance RMS' );
bar( [ batch.sigma_rms , seq.sigma_rms ] );
set( gca , 'XTickLabel' , labels );
ylabel( 'position \sigma_{rms} [m]' );
legend( 'Batch' , 'Sequential' , 'Location' , 'best' );
grid on
title( 'Position Covariance RMS vs Station P0 Variance' );

figure( 'Name' , 'P0 Sensitivity: Position Covariance Trace' );
bar( [ batch.trace , seq.trace ] );
set( gca , 'XTickLabel' , labels );
ylabel( 'trace(P_{pos}) [m^2]' );
legend( 'Batch' , 'Sequential' , 'Location' , 'best' );
grid on
title( 'Position Covariance Trace vs Station P0 Variance' );

save project_sensitivity_p0_results stationVars batch seq

fprintf( '\nSaved results to project_sensitivity_p0_results.mat\n' );
