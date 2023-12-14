function cs = cos_sim(x, y)
% Cosine similarity between two vectors x and y. Author: Jonathan Chien.

x_ = x ./ vecnorm(x);
y_ = y ./ vecnorm(y);

cs = dot(x_, y_);

end