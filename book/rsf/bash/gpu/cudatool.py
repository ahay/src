def CUDANVCCStaticObjectEmitter(target, source, env):
    import os
    import SCons.Defaults
    tgt, src = SCons.Defaults.StaticObjectEmitter(target, source, env)
    for file in src:
        linkfn = os.path.basename(
            os.path.splitext(src[0].rstr())[0]) + '.linkinfo'
        env.Clean(src, linkfn)
    return tgt, src
def CUDANVCCSharedObjectEmitter(target, source, env):
    import os
    import SCons.Defaults
    tgt, src = SCons.Defaults.SharedObjectEmitter(target, source, env)
    for file in src:
        linkfn = os.path.basename(
            os.path.splitext(src[0].rstr())[0]) + '.linkinfo'
        env.Clean(src, linkfn)
    return tgt, src

def find_paths(paths):
    """
    @param paths: paths to geuss.
    @type paths: list
    @return: the path found or None (not found).
    @rtype: str
    """
    import os
    for path in paths:
        if os.path.isdir(path):
            return path
    return ''

def determine_paths(env):
    """
    Fill the 'CUDA_TOOLKIT_PATH' and 'CUDA_SDK_PATH' into environment if they
    were not there.

    @return: the paths.
    @rtype: tuple
    """
    import sys
    import os
    from warnings import warn
    home = os.environ.get('HOME', '')
    programfiles = os.environ.get('PROGRAMFILES', '')
    homedrive = os.environ.get('HOMEDRIVE', '')

    # find CUDA Toolkit path and set CUDA_TOOLKIT_PATH.
    cudaToolkitPath = os.environ.get('CUDA_TOOLKIT_PATH', '')
    if not cudaToolkitPath:
        paths = [
            '/'.join([home, 'NVIDIA_CUDA_TOOLKIT']),
            '/'.join([home, 'Apps', 'NVIDIA_CUDA_TOOLKIT']),
            '/'.join([home, 'Apps', 'CudaToolkit']),
            '/'.join([home, 'Apps', 'CudaTK']),
            '/'.join(['/usr', 'local', 'NVIDIA_CUDA_TOOLKIT']),
            '/'.join(['/usr', 'local', 'CUDA_TOOLKIT']),
            '/'.join(['/usr', 'local', 'cuda_toolkit']),
            '/'.join(['/usr', 'local', 'CUDA']),
            '/'.join(['/usr', 'local', 'cuda']),
            '/'.join(['/Developer', 'NVIDIA CUDA TOOLKIT']),
            '/'.join(['/Developer', 'CUDA TOOLKIT']),
            '/'.join(['/Developer', 'CUDA']),
            '/'.join([programfiles, 'NVIDIA Corporation',
                'NVIDIA CUDA TOOLKIT']),
            '/'.join([programfiles, 'NVIDIA Corporation', 'NVIDIA CUDA']),
            '/'.join([programfiles, 'NVIDIA Corporation', 'CUDA TOOLKIT']),
            '/'.join([programfiles, 'NVIDIA Corporation', 'CUDA']),
            '/'.join([programfiles, 'NVIDIA', 'NVIDIA CUDA TOOLKIT']),
            '/'.join([programfiles, 'NVIDIA', 'NVIDIA CUDA']),
            '/'.join([programfiles, 'NVIDIA', 'CUDA TOOLKIT']),
            '/'.join([programfiles, 'NVIDIA', 'CUDA']),
            '/'.join([programfiles, 'CUDA TOOLKIT']),
            '/'.join([programfiles, 'CUDA']),
            '/'.join([homedrive, 'CUDA TOOLKIT']),
            '/'.join([homedrive, 'CUDA']),
        ]
        cudaToolkitPath = find_paths(paths)
        if cudaToolkitPath:
            sys.stdout.write(
                'scons: CUDA Toolkit found in %s\n' % cudaToolkitPath)
        #else:
        #    warn('Cannot find the CUDA Toolkit path. '
        #        'Please set it to CUDA_TOOLKIT_PATH environment variable.')
    env['CUDA_TOOLKIT_PATH'] = cudaToolkitPath

    # find CUDA SDK path and set CUDA_SDK_PATH.
    cudaSDKPath = os.environ.get('CUDA_SDK_PATH', '')
    if not cudaSDKPath:
        paths = [
            '/'.join([home, 'NVIDIA_GPU_Computing_SDK', 'C']),
            '/'.join([home, 'NVIDIA_GPU_Computing_SDK']),
            '/'.join([home, 'NVIDIA_CUDA_SDK', 'C']),
            '/'.join([home, 'NVIDIA_CUDA_SDK']),
            '/'.join([home, 'opt', 'cudasdk', 'C']),
            '/'.join([home, 'opt', 'NVIDIA_GPU_Computing_SDK', 'C']),
            '/'.join([home, 'opt', 'NVIDIA_CUDA_SDK', 'C']),
            '/'.join([home, 'Apps', 'NVIDIA_CUDA_SDK']),
            '/'.join([home, 'Apps', 'CudaSDK']),
            '/'.join(['/usr', 'local', 'NVIDIA_CUDA_SDK']),
            '/'.join(['/usr', 'local', 'CUDASDK']),
            '/'.join(['/usr', 'local', 'cuda_sdk']),
            '/'.join(['/Developer', 'NVIDIA CUDA SDK']),
            '/'.join(['/Developer', 'CUDA SDK']),
            '/'.join(['/Developer', 'CUDA']),
            '/'.join([programfiles, 'NVIDIA Corporation', 'NVIDIA CUDA SDK']),
            '/'.join([programfiles, 'NVIDIA', 'NVIDIA CUDA SDK']),
            '/'.join([programfiles, 'NVIDIA CUDA SDK']),
            '/'.join([programfiles, 'CudaSDK']),
            '/'.join([homedrive, 'NVIDIA CUDA SDK']),
            '/'.join([homedrive, 'CUDA SDK']),
            '/'.join([homedrive, 'CUDA', 'SDK']),
        ]
        cudaSDKPath = find_paths(paths)
        if cudaSDKPath:
            sys.stdout.write(
                'scons: CUDA SDK found in %s \n' % cudaSDKPath)
        #else:
        #    warn('Cannot find the CUDA SDK path. '
        #        'Please set it to CUDA_SDK_PATH environment variable.')
    env['CUDA_SDK_PATH'] = cudaSDKPath

    # return.
    return cudaToolkitPath, cudaSDKPath

def generate(env):
    """
    In order to use this tool, user must have CUDA_SDK_PATH and
    CUDA_TOOLKIT_PATH environmental variables defined.
    """
    import os
    import SCons.Tool
    import SCons.Scanner.C
    cudaToolkitPath, cudaSDKPath = determine_paths(env)

    # scanners and builders.
    CUDAScanner = SCons.Scanner.C.CScanner()
    staticObjBuilder, sharedObjBuilder = SCons.Tool.createObjBuilders(env);
    staticObjBuilder.add_action('.cu', '$STATICNVCCCMD')
    staticObjBuilder.add_emitter('.cu', CUDANVCCStaticObjectEmitter)
    sharedObjBuilder.add_action('.cu', '$SHAREDNVCCCMD')
    sharedObjBuilder.add_emitter('.cu', CUDANVCCSharedObjectEmitter)
    SCons.Tool.SourceFileScanner.add_scanner('.cu', CUDAScanner)

    # build commands.
    env['STATICNVCCCMD'] = ' '.join([
        '$NVCC',
        '$NVCCINC',
        '$NVCCFLAGS',
        '$STATICNVCCFLAGS',
        '-c $SOURCES',
    ])
    env['SHAREDNVCCCMD'] = ' '.join([
        '$NVCC',
        '$NVCCINC',
        '$NVCCFLAGS',
        '$SHAREDNVCCFLAGS',
        '$ENABLESHAREDNVCCFLAG',
        '-c $SOURCES',
    ])

    # compiler.
    env['NVCC'] = 'nvcc'
    env.PrependENVPath('PATH', '/'.join([cudaToolkitPath, 'bin']))

    # includes.
    env.Append(NVCCINC=' '.join(['-I',
        '/'.join([cudaSDKPath, 'common', 'inc']),
    ]))
    env.Append(CPPPATH=[
        '/'.join([cudaSDKPath, 'common', 'inc']),
        '/'.join([cudaToolkitPath, 'include']),
    ])

    # NVCC compiler flags.
    env['NVCCFLAGS'] = ''
    env['STATICNVCCFLAGS'] = ''
    env['SHAREDNVCCFLAGS'] = ''
    env['ENABLESHAREDNVCCFLAG'] = '-shared'

    # libraries.
    if env['PLATFORM'] == 'posix':
        cudaSDKSubLibDir = 'linux'
    elif env['PLATFORM'] == 'darwin':
        cudaSDKSubLibDir = 'darwin'
    else:
        cudaSDKSubLibDir = ''
    env.Append(LIBPATH=[
        '/'.join([cudaSDKPath, 'lib']),
        '/'.join([cudaSDKPath, 'common', 'lib', cudaSDKSubLibDir]),
        '/'.join([cudaToolkitPath, 'lib64']),
        '/'.join([cudaToolkitPath, 'lib']),
    ])

def exists(env):
    return env.Detect('nvcc')

