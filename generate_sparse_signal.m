function signal_matrix = generate_sparse_signal(signal_length, num_measurements, sparsity_type, sparsity_level, num_blocks, block_size, isolated_ratio)
    % Generate a sparse MMV signal with different types of sparsity in MATLAB.
    %
    % Parameters:
    % signal_length: int, total length of each measurement vector (number of rows).
    % num_measurements: int, number of measurement vectors (number of columns).
    % sparsity_type: string, type of sparsity ('random', 'block', or 'hybrid').
    % sparsity_level: float, fraction of non-zero rows for 'random' sparsity.
    % num_blocks: int, number of non-zero blocks for 'block' sparsity.
    % block_size: int, size of each non-zero block for 'block' and 'hybrid' sparsity.
    % isolated_ratio: float, fraction of isolated non-zero rows for 'hybrid' sparsity.
    %
    % Returns:
    % signal_matrix: matrix of size (signal_length x num_measurements) representing the sparse MMV signal.

    % Initialize a zero matrix of given dimensions
    signal_matrix = zeros(signal_length, num_measurements);
    
    if strcmp(sparsity_type, 'Random')
        % Randomly choose rows to be non-zero across all measurement vectors
        num_nonzeros = round(sparsity_level * signal_length);
        non_zero_rows = randperm(signal_length, num_nonzeros);
        signal_matrix(non_zero_rows, :) = sqrt(.5)*(randn(num_nonzeros, num_measurements)+1j*randn(num_nonzeros, num_measurements));
        
    elseif strcmp(sparsity_type, 'Block')
        % Generate non-zero blocks in rows that are shared across all measurement vectors
        block_positions = randperm(signal_length - block_size + 1, num_blocks);
        for i = 1:num_blocks
            pos = block_positions(i);
            signal_matrix(pos:pos + block_size - 1, :) = sqrt(.5)*(randn(block_size, num_measurements)+1j*randn(block_size, num_measurements));
        end
        
    elseif strcmp(sparsity_type, 'Hybrid')
        % Hybrid sparsity with both isolated and block non-zero rows
        
        % Determine number of isolated non-zero rows
        num_isolated = round(isolated_ratio * sparsity_level * signal_length);
        isolated_rows = randperm(signal_length, num_isolated);
        signal_matrix(isolated_rows, :) = sqrt(.5)*(randn(num_isolated, num_measurements)+1j*randn(num_isolated, num_measurements));
        
        % Determine number of blocks and create blocks
        num_blocks = round((1 - isolated_ratio) * sparsity_level * signal_length / block_size);
        block_positions = randperm(signal_length - block_size + 1, num_blocks);
        for i = 1:num_blocks
            pos = block_positions(i);
            signal_matrix(pos:pos + block_size - 1, :) = sqrt(.5)*(randn(block_size, num_measurements)+1j*randn(block_size, num_measurements));
        end
        
    else
        error("Invalid sparsity_type. Choose from 'Random', 'Block', or 'Hybrid'.")
    end
end

