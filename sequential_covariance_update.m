%-------------------------------------------------------------------------%
%    sequential_covariance_update
%
%  DESCRIPTION:
%     Measurement-update for the sequential covariance matrix.
%
%     The textbook conventional form P = (I-KH)Pbar is algebraically
%     correct but often loses positive semi-definiteness in finite
%     precision, which later produces complex ellipsoid axes.  When
%     S.useJoseph is true, the Joseph form is used instead.
%
%  MAE 182 -- Spacecraft Guidance and Navigation
%-------------------------------------------------------------------------%
function P_k = sequential_covariance_update( S , P_bar , K , Htilde , I18 )

if S.useJoseph
    IKH = I18 - K * Htilde;
    P_k = IKH * P_bar * IKH' + K * S.R * K';
else
    Syy = Htilde * P_bar * Htilde' + S.R;
    P_k = P_bar - K * Syy * K';
end

if S.symmetrizeCov
    P_k = 0.5 * ( P_k + P_k' );
end

P_k = real( P_k );

end
